#!/usr/bin/python3
from sra_fastq_dump_qt import *
from qt5_gui_sra import Graph_UI
from PyQt5.QtWidgets import QApplication
from Bio import Entrez
import GEOparse
import sys



def main():
    app = QApplication(sys.argv)
    qt5_data = Graph_UI()
    app.exec_()
    gs_text, series, s3_text, output_path, email = read_qt5_input(qt5_data)
    py3 = sys.version_info[0] > 2 #creates boolean value for test that Python major version > 2

    try:
        os.mkdir(output_path)
        response = 'N'
    except OSError:
        if py3:
          response = input("Directory already exists. Do you still want to continue (Y/N): ")
        else:
          response = raw_input("Directory already exists. Do you still want to continue (Y/N): ")
        if response != 'Y':
            sys.exit("Canceled, please rerun with different Directory or path.")
    if gs_text != '':
        sra_series_fetch(gs_text, series, s3_text, output_path, email, response)
if __name__ == '__main__':
    main()
