from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys

def get_parser():
    """Parser command line args."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-geo", "--sra",
                        dest="sra",
                        type=str,
                        required =True,
                        help="GSM or GSE geo number")
    parser.add_argument("-s3", "--s3_bucket",
                        dest="s3",
                        type=str,
                        default='',
                        help="s3 bucket if desired.")
    parser.add_argument("-out", "-o",
                        dest="directoryLabel",
                        type=str,
                        default='',
                        help="path to download files.")
    parser.add_argument("-e", "-email",
                        dest="email",
                        type=str,
                        default='',
                        help="Email for ncbi tools.")
    parser.add_argument("-local",
                        dest="local_manifest",
                        default=False,
                        action='store_true',
                        help="Make manifest of files for GSE.")
    parser.add_argument("-s3_manifest",
                        dest="s3_manifest",
                        default=False,
                        action='store_true',
                        help="Make manifest of files for GSE.")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        from qt5_gui_sra import Graph_UI
        from PyQt5.QtWidgets import QApplication
        app = QApplication(sys.argv)
        qt5_data = Graph_UI()
        app.exec_()
        return ('qt', qt5_data)
    args = parser.parse_args()
    return ('arg',args)
