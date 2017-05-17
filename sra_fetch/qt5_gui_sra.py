import sys
from PyQt5.QtWidgets import (QWidget,QMainWindow, QTextEdit,
    QAction, QFileDialog, QApplication, QPushButton, QLineEdit, QGridLayout, QCheckBox, QInputDialog, QLabel)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt, QCoreApplication


class Graph_UI(QWidget):

    def __init__(self):
        super().__init__()

        self.initUI()


    def initUI(self):
        #setup grid layout
        grid = QGridLayout()
        grid.setSpacing(10)

        #get SRA number
        self.SRA_label = QLabel("SRA GSM or GSE #: ")
        self.sra = QLineEdit(self)
        self.sra.setFixedWidth(200)
        self.sra.setText('')

        #enter s3 path if upload is desired
        self.s3_label = QLabel('Upload to S3 path:')
        self.s3 = QLineEdit(self)
        self.s3.setFixedWidth(400)
        self.s3.setText('')

        #choose directory for output if not current directory
        self.btn = QPushButton('Select Directory for Download', self)
        self.btn.clicked.connect(self.fileDialog)
        self.dir_select_label = QLabel("Select the directory for download. ")
        self.directoryLabel = QLineEdit(self)
        self.directoryLabel.setFixedWidth(300)
        self.directoryLabel.setText('')

        #enter email for NCBI Entrez
        self.email_label = QLabel('Email for NCBI Entrez:')
        self.email = QLineEdit(self)
        self.email.setFixedWidth(200)
        self.email.setText('')

        self.qbtn = QPushButton('Run', self)
        self.qbtn.clicked.connect(self.get_sra)
        self.qbtn.clicked.connect(self.get_s3)
        self.qbtn.clicked.connect(self.get_email)
        self.qbtn.clicked.connect(QCoreApplication.instance().quit)
        self.qbtn.resize(self.qbtn.sizeHint())
        self.qbtn.move(50, 50)



        grid.addWidget(self.SRA_label, 0,0)
        grid.addWidget(self.sra, 1,0)
        grid.addWidget(self.dir_select_label, 0,1)
        grid.addWidget(self.directoryLabel, 1,1)
        grid.addWidget(self.btn, 2,1)
        grid.addWidget(self.s3_label, 3,0)
        grid.addWidget(self.s3, 4,0)
        grid.addWidget(self.email_label, 3,1)
        grid.addWidget(self.email, 4,1)
        grid.addWidget(self.qbtn, 5,0)

        self.setLayout(grid)

        #self.setGeometry(300, 300, 350, 300)
        self.setWindowTitle('SRA Downloader')
        self.show()


    def fileDialog(self):
        options = QFileDialog.DontResolveSymlinks | QFileDialog.ShowDirsOnly
        directory = QFileDialog.getExistingDirectory(self,
                "QFileDialog.getExistingDirectory()",
                self.directoryLabel.text(), options=options)
        if directory:
            self.directoryLabel.setText(directory)

    def get_sra(self):
        text = self.sra.text()
        self.sra.setText(text)
        print('SRA input: ',text)

    def get_s3(self):
        s3text = self.s3.text()
        self.s3.setText(s3text)
        if s3text != '':
            print('s3 upload path: ', s3text)

    def get_email(self):
        email = self.email.text()
        self.email.setText(email)
        if email != '':
            print('email for NCBI: ', email)
