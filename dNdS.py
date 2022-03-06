# importing required modules
import os
import sys
import math
import pandas as pd
import tkinter as tk
import numpy as np
import colorama
import requests
from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog, messagebox
from pandastable import Table, TableModel
from PIL import ImageTk, Image
from io import BytesIO

from values import Values

# to render colored text
os.system("")  # enables ansi escape characters

COLOR = {
    "HEADER": "\033[95m",
    "BLUE": "\033[94m",
    "GREEN": "\033[92m",
    "RED": "\033[91m",
    "ENDC": "\033[0m",
    "BOLD": "\033[1m",
}


class Page(tk.Toplevel):
    def __init__(self, *args, **kwargs):
        super(Page, self).__init__(*args, **kwargs)
        # self.geometry('1500x500')
        self.title("Calculate dN/dS value")
        self.grab_set()

        # labels
        self.step1_opt1_lbl = tk.Label(
            self,
            text="Step I: Upload gene sequences of 1st organism in FASTA format",
            font=("Arial", 12, "bold"),
            fg="brown",
        )
        self.step1_opt2_lbl = tk.Label(
            self,
            text="Paste gene sequences of 1st organism in FASTA format",
            font=("Arial", 12, "bold"),
            fg="blue",
        )
        self.step2_opt1_lbl = tk.Label(
            self,
            text="Step 2: Upload gene sequences of 2nd organism in FASTA format",
            font=("Arial", 12, "bold"),
            fg="brown",
        )
        self.step2_opt2_lbl = tk.Label(
            self,
            text="Paste gene sequences of 2nd organism in FASTA format",
            font=("Arial", 12, "bold"),
            fg="blue",
        )

        # process button
        style = Style()
        style.configure(
            "TButton", font=("arial", 12, "bold"), borderwidth="4", relief=RAISED
        )
        style.map(
            "TButton",
            foreground=[("active", "!disabled", "green")],
            background=[("active", "black")],
        )

        self.process = Button(
            self, text="Process", style="TButton", command=self.process
        )

        # or labels
        self.or1_lbl = tk.Label(self, text="Or", font=("Arial", 11, "bold"))
        self.or2_lbl = tk.Label(self, text="Or", font=("Arial", 11, "bold"))

        # upload buttons
        self.seq1_btn = Button(
            self, text="Upload", style="TButton", command=self.get_file1
        )
        self.seq2_btn = Button(
            self, text="Upload", style="TButton", command=self.get_file2
        )
        self.uploaded_file1_lbl = tk.Label(self, text="")
        self.uploaded_file2_lbl = tk.Label(self, text="")

        # textboxes
        self.seq1_txtbox = tk.Text(self, height=10, width=45)
        self.seq2_txtbox = tk.Text(self, height=10, width=45)

        # grid packing
        self.step1_opt1_lbl.grid(row=0, column=0)
        self.seq1_btn.grid(row=1, column=0)
        self.uploaded_file1_lbl.grid(row=2, column=0)
        self.or1_lbl.grid(row=3, column=0)
        self.step1_opt2_lbl.grid(row=4, column=0)
        self.seq1_txtbox.grid(row=5, column=0, padx=10, pady=10)

        self.step2_opt1_lbl.grid(row=0, column=1, padx=25)
        self.seq2_btn.grid(row=1, column=1)
        self.uploaded_file2_lbl.grid(row=2, column=1)
        self.or2_lbl.grid(row=3, column=1)
        self.step2_opt2_lbl.grid(row=4, column=1)
        self.seq2_txtbox.grid(row=5, column=1, padx=10, pady=10)

        self.process.grid(row=6, column=0, columnspan=2)

    def get_file1(self):
        self.file1_path = filedialog.askopenfile(
            mode="r", filetypes=[("Text Files", ".txt")]
        )
        self.uploaded_file1_lbl["text"] = "Selected file: " + self.file1_path.name

    def get_file2(self):
        self.file2_path = filedialog.askopenfile(
            mode="r", filetypes=[("Text Files", ".txt")]
        )
        self.uploaded_file2_lbl["text"] = "Selected file: " + self.file2_path.name

    def process(self):
        """Process content from files uploaded, if files not uploaded, process from textbox"""
        try:
            content1 = [line for line in self.file1_path.read().splitlines() if line]
            content2 = [line for line in self.file2_path.read().splitlines() if line]
        except:
            content1 = [
                line for line in self.seq1_txtbox.get("1.0", "end").split("\n") if line
            ]
            content2 = [
                line for line in self.seq2_txtbox.get("2.0", "end").split("\n") if line
            ]

        codon = Codon(content1, content2)
        self.display_df(codon.main_df)

    def display_df(self, df):
        """Display the final result dataframe that is formed"""
        # create child window
        win = Toplevel()
        win.geometry("900x800")
        win.title("Details per codon sequence")
        win.grab_set()
        # using pandastable for better formatting
        pt = Table(
            win, dataframe=df.fillna("Math error"), showtoolbar=True, showstatusbar=True
        )
        pt.show()


class Codon:
    def __init__(self, file1, file2):
        values = Values()
        self.AA3 = values.AA3
        self.AA1 = values.AA1
        self.Codon = values.Codon
        self.SynSite = values.SynSite
        self.NonSynSite = values.NonSynSite

        """ Create dictionary maps of the above information
        """
        self.codon_amino_map = {key: val for key, val in zip(self.Codon, self.AA1)}
        self.codon_synsite_map = {
            key: val for key, val in zip(self.Codon, self.SynSite)
        }
        self.codon_nonsynsite_map = {
            key: val for key, val in zip(self.Codon, self.NonSynSite)
        }

        self.file1 = file1
        self.file2 = file2

        self.out_file_name = "output.csv"

        self.read_files()
        self.calculate_syn_nonsyn_sites()
        self.calculate_syn_nonsyn_changes()
        self.calculate_pN_pS()
        self.calculate_dN_dS()

        self.write_values_to_file()


    def read_files(self):
        """Read files to make a list of genes and their codon sequences for each organism"""
        self.info1 = [line for line in self.file1 if line.startswith(">")]
        self.refSeq1 = [line for line in self.file1 if not line.startswith(">")]

        self.info2 = [line for line in self.file2 if line.startswith(">")]
        self.refSeq2 = [line for line in self.file2 if not line.startswith(">")]

        # check for codon length mismatch - throw error
        for i, (seq1, seq2) in enumerate(zip(self.refSeq1, self.refSeq2)):
            if len(seq1) != len(seq2):
                messagebox.showerror(
                    "Error",
                    f"Length of {i+1}th sequences are not equal!\nUsing min length for current calculation!",
                )

    def calculate_syn_nonsyn_sites(self):
        """calculate synonymous and non-synonymous SITES per codon"""
        self.total_synsites_per_codon = []
        self.total_nonsynsites_per_codon = []

        for k in range(len(self.refSeq1)):
            synsites1 = synsites2 = nonsynsites1 = nonsynsites2 = 0
            min_length = min(len(self.refSeq1[k]), len(self.refSeq2[k]))
            for i in range(0, min_length, 3):
                code1 = self.refSeq1[k][i : i + 3].upper().replace("T", "U")
                code2 = self.refSeq2[k][i : i + 3].upper().replace("T", "U")

                synsites1 += self.codon_synsite_map[code1]
                synsites2 += self.codon_synsite_map[code2]
                nonsynsites1 += self.codon_nonsynsite_map[code1]
                nonsynsites2 += self.codon_nonsynsite_map[code2]

            self.total_synsites_per_codon.append((synsites1 + synsites2) / 2)
            self.total_nonsynsites_per_codon.append((nonsynsites1 + nonsynsites2) / 2)

    def calculate_syn_nonsyn_changes(self):
        """calculate synonymous and non-synonymous CHANGES per codon"""
        self.total_synmut_per_codon = []
        self.total_nonsynmut_per_codon = []

        for k in range(len(self.refSeq1)):
            total_syn_mut = total_nonsyn_mut = 0
            min_length = min(len(self.refSeq1[k]), len(self.refSeq2[k]))
            for i in range(0, min_length, 3):
                code1 = self.refSeq1[k][i : i + 3].upper().replace("T", "U")
                code2 = self.refSeq2[k][i : i + 3].upper().replace("T", "U")

                if code1 != code2:
                    synMut = nonSynMut = total = 0
                    self.paths = []
                    self.find_paths(code1, [code1], code1, code2)

                    for path in self.paths:
                        synMut_per_path = nonSynMut_per_path = 0

                        # if first codon is STOP codon, ignore path
                        if self.codon_amino_map[path[0]] == "X":
                            continue
                        for ind in range(len(path) - 1):
                            # if codon is STOP codon, ignore path
                            if self.codon_amino_map[path[ind + 1]] == "X":
                                synMut_per_path = nonsynMut_per_path = 0
                                break
                            elif (
                                self.codon_amino_map[path[ind]]
                                == self.codon_amino_map[path[ind + 1]]
                            ):
                                synMut_per_path += 1
                            else:
                                nonSynMut_per_path += 1

                        synMut += synMut_per_path
                        nonSynMut += nonSynMut_per_path

                    total = synMut + nonSynMut

                    # when first codon triplet for the only path is a STOP codon
                    if total > 0:
                        total_syn_mut += (synMut / total) * (len(self.paths[0]) - 1)
                        total_nonsyn_mut += (nonSynMut / total) * (
                            len(self.paths[0]) - 1
                        )

            self.total_synmut_per_codon.append(total_syn_mut)
            self.total_nonsynmut_per_codon.append(total_nonsyn_mut)

    def find_paths(self, str1, str2, main_str, out):
        """A recursive function to find all paths from a given codon to another"""
        if str1 == out:
            self.paths.append(str2)
            return
        for i in range(len(main_str)):
            if str1[i] != out[i]:
                new_str = list(str1)
                new_str[i] = out[i]
                new_str = "".join(new_str)
                self.find_paths(new_str, str2 + [new_str], main_str, out)

    def calculate_pN_pS(self):
        self.pS_per_codon = [
            x / y
            for x, y in zip(self.total_synmut_per_codon, self.total_synsites_per_codon)
        ]
        self.pN_per_codon = [
            x / y
            for x, y in zip(
                self.total_nonsynmut_per_codon, self.total_nonsynsites_per_codon
            )
        ]

        try:
            self.pNpS = [i / j for i, j in zip(self.pN_per_codon, self.pS_per_codon)]
        except:
            messagebox.showerror(
                "Error",
                "Error in calculating pN/pS value due to ZeroDivision Error! Please recheck your codon sequences.",
            )

    def calculate_dN_dS(self):
        def util_func(i):
            """The codons for which the calculation does not render valid values are made NaN"""
            try:
                return -(3 / 4) * math.log(1 - ((4 / 3) * i))
            except:
                return np.nan

        dS_per_codon = [util_func(i) for i in self.pS_per_codon]
        dN_per_codon = [util_func(i) for i in self.pN_per_codon]

        # print("dN_per_codon", dN_per_codon)
        # print("dS_per_codon", dS_per_codon)

        self.dNdS = [i / j for i, j in zip(dN_per_codon, dS_per_codon)]

    def write_values_to_file(self):
        self.main_df = pd.DataFrame(
            {
                "Gene": self.info1,
                "Size (codon)": [len(codon) for codon in self.refSeq1],
                "No. Syn Mutation": self.total_synmut_per_codon,
                "No. NonSyn Mutation": self.total_nonsynmut_per_codon,
                "SynSites": self.total_synsites_per_codon,
                "NonSynSites": self.total_nonsynsites_per_codon,
                "pN/pS": self.pNpS,
                "dN/dS": self.dNdS,
                "Remark": [
                    ""
                    if len(c1) == len(c2)
                    else "Length mismatch. Considering length of first sequence"
                    for c1, c2 in zip(self.refSeq1, self.refSeq2)
                ],
            }
        )

        self.main_df.to_csv(self.out_file_name, index=False, na_rep="NaN")


class AboutPage(tk.Toplevel):
    def __init__(self, *args, **kwargs):
        super(AboutPage, self).__init__(*args, **kwargs)
        self.title("About the authors")
        self.geometry("450x200")
        # labels
        self.text = tk.Label(
            self,
            text="CUBCal is a software being developed under the guidance of Dr.S.S Sathapathy, \nan associate professor of computer science and engeneering, Tezpur University.\nThis tool will provide the best alternative to the already existing Codonw software,\nwidely used for codon usage analysis and still under contruction",
        )
        self.text.grid(row=0, column=0)


if __name__ == "__main__":
    # root page
    root = Tk()
    root.title("Main Menu")
    root.geometry("1000x400")

    # insert TU logo 
    image_url = "https://bihog.com/wp-content/uploads/2020/01/Tezpur-University-1200x814.png"
    response = requests.get(image_url)
    img = PhotoImage(data=response.content)

    # Create a Canvas
    canvas = Canvas(root, width=1000, height=400)
    canvas.pack(fill=BOTH, expand=True)

    # Add Image inside the Canvas
    canvas.create_image(0, 0, image=img, anchor="nw")

    # Function to resize the window
    def resize_image(e):
        global image, resized, image2
        # open image to resize it
        image = Image.open(BytesIO(response.content))
        # resize the image with width and height of root
        resized = image.resize((e.width, e.height), Image.ANTIALIAS)

        image2 = ImageTk.PhotoImage(resized)
        canvas.create_image(0, 0, image=image2, anchor="nw")

    # Bind the function to configure the parent window
    root.bind("<Configure>", resize_image)

    # content = 'CUBCal is a software being developed under the guidance of Dr.S.S Sathapathy, \nan associate professor of computer science and engeneering, Tezpur University.\nThis tool will provide the best alternative to the already existing Codonw software,\nwidely used for codon usage analysis and still under contruction'
    # text = Text(root)
    # text.insert(INSERT, content)
    # text.pack(expand= 1, fill= BOTH)
    # text.configure(font=("Times", 16, "bold"))
    # text.tag_config("start", foreground="red")
    # text.tag_add("start", "1.0", "1.6")
    # text.tag_add("start", "1.59", "1.76")

    # Creating Menubar
    menubar = Menu(root)

    # Adding File Menu and commands
    file = Menu(menubar, tearoff=0)
    menubar.add_cascade(label="File", menu=file)
    file.add_command(label="New File", command=None)
    file.add_command(label="Open...", command=None)
    file.add_command(label="Save", command=None)
    file.add_separator()
    file.add_command(label="Exit", command=root.destroy)

    # Adding Edit Menu and commands
    edit = Menu(menubar, tearoff=0)
    menubar.add_cascade(label="Edit", menu=edit)
    edit.add_command(label="Cut", command=None)
    edit.add_command(label="Copy", command=None)
    edit.add_command(label="Paste", command=None)
    edit.add_command(label="Select All", command=None)
    edit.add_separator()
    file.add_command(label="Exit", command=root.destroy)

    # Adding Help Menu
    page = lambda: Page()
    help_ = Menu(menubar, tearoff=0)
    menubar.add_cascade(label="Selection", menu=help_)
    help_.add_command(label="dN/dS", command=page)
    help_.add_command(label="Demo", command=None)
    help_.add_separator()
    help_.add_command(label="About Tk", command=None)

    # Codon usage index
    cui = Menu(menubar, tearoff=0)
    menubar.add_cascade(label="Codon Usage Index", menu=cui)
    cui.add_command(label="Codon Adaptation Index", command=None)
    # cui.add_command(label ='Copy', command = None)
    # cui.add_command(label ='Paste', command = None)
    # cui.add_command(label ='Select All', command = None)

    # about
    about_page = lambda: AboutPage()
    about = Menu(menubar, tearoff=0)
    menubar.add_cascade(label="About", menu=about)
    about.add_command(label="About the creators", command=about_page)

    # starting the page
    root.config(menu=menubar)
    root.mainloop()
