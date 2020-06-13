#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This is AGEseqPy
"""

import sys
from sys import argv
import glob, string
import os
import os.path
import re
from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_AGEseq(object):

    def setupUi(self, AGEseq):
        AGEseq.setObjectName("AGEseq")
        AGEseq.resize(540, 386)
        self.centralwidget = QtWidgets.QWidget(AGEseq)
        self.centralwidget.setObjectName("centralwidget")
        self.file_button = QtWidgets.QPushButton(self.centralwidget)
        self.file_button.setGeometry(QtCore.QRect(50, 40, 151, 41))
        font = QtGui.QFont()
        font.setFamily("Agency FB")
        font.setPointSize(9)
        self.file_button.setFont(font)
        self.file_button.setObjectName("file_button")
        self.goalfile_button = QtWidgets.QPushButton(self.centralwidget)
        self.goalfile_button.setGeometry(QtCore.QRect(50, 110, 151, 41))
        font = QtGui.QFont()
        font.setFamily("Agency FB")
        font.setPointSize(9)
        self.goalfile_button.setFont(font)
        self.goalfile_button.setObjectName("goalfile_button")
        self.run_button = QtWidgets.QPushButton(self.centralwidget)
        self.run_button.setGeometry(QtCore.QRect(50, 270, 151, 41))
        font = QtGui.QFont()
        font.setFamily("Agency FB")
        font.setPointSize(9)
        self.run_button.setFont(font)
        self.run_button.setObjectName("run_button")
        self.file_input = QtWidgets.QTextBrowser(self.centralwidget)
        self.file_input.setGeometry(QtCore.QRect(230, 40, 281, 41))
        self.file_input.setObjectName("file_input")
        self.goalfile_input = QtWidgets.QTextBrowser(self.centralwidget)
        self.goalfile_input.setGeometry(QtCore.QRect(230, 110, 281, 41))
        self.goalfile_input.setObjectName("goalfile_input")
        self.run_output = QtWidgets.QTextBrowser(self.centralwidget)
        self.run_output.setGeometry(QtCore.QRect(230, 250, 281, 81))
        self.run_output.setObjectName("run_output")
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(50, 180, 151, 41))
        font = QtGui.QFont()
        font.setFamily("Agency FB")
        font.setPointSize(9)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.textBrowser = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser.setGeometry(QtCore.QRect(230, 180, 281, 41))
        self.textBrowser.setObjectName("textBrowser")
        AGEseq.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(AGEseq)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 540, 18))
        self.menubar.setObjectName("menubar")
        AGEseq.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(AGEseq)
        self.statusbar.setObjectName("statusbar")
        AGEseq.setStatusBar(self.statusbar)

        self.retranslateUi(AGEseq)
        QtCore.QMetaObject.connectSlotsByName(AGEseq)

    def retranslateUi(self, AGEseq):
        _translate = QtCore.QCoreApplication.translate
        AGEseq.setWindowTitle(_translate("AGEseq", "AGEseq"))
        self.file_button.setText(_translate("AGEseq", "选择文件夹"))
        self.goalfile_button.setText(_translate("AGEseq", "选择目标文件"))
        self.run_button.setText(_translate("AGEseq", "运行"))
        self.pushButton.setText(_translate("AGEseq", "选择输出文件"))



class MainUi(QtWidgets.QMainWindow, Ui_AGEseq):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_AGEseq.__init__(self)
        self.setupUi(self)
        self.file_button.clicked.connect(self.selfile)
        self.goalfile_button.clicked.connect(self.selgoalfile)
        self.pushButton.clicked.connect(self.selfileout)
        self.run_button.clicked.connect(self.run)
        self.fileDialog = QtWidgets.QFileDialog(self)

    def selfile(self):
        self.selfile_input = self.fileDialog.getExistingDirectory()
        self.file_input.setText("".join(self.selfile_input)+"\n")
        print("input"+"".join(self.selfile_input))

    def selgoalfile(self):
        self.selgoalfile_input,_filter = self.fileDialog.getOpenFileName()
        self.goalfile_input.setText("".join(self.selgoalfile_input))
        print("input"+"".join(self.selgoalfile_input))

    def selfileout(self):
        self.file_out, _filter = self.fileDialog.getSaveFileName()
        self.textBrowser.setText("".join(self.file_out))
        print("output"+"".join(self.file_out))

    def run(self):
        out = "blat Job Starts"+"\n"
        self.run_output.setText(out)
        script, dat_dir, file_design = argv, self.selfile_input, self.selgoalfile_input

        blat_dir = ''  # working directory or PATH

        ##########################
        # setting can be changed here

        # setting for reports
        mismatch_cutoff = 0.1  # mismatch rate to filter low quality alignment, default 0.1 (10 %)
        min_cutoff = 0  # cutoff to filter reads with low abundance, default  0
        wt_like_report = 20  # report top xx WT like records, default 20
        indel_report = 50  # report top xx records with indel, default  50
        remove_files = 1  # keep (0) or delete (1) intermediate files, default = 1

        remove = 'del'

        ################
        # input

        design_fas = file_design + '_DESIGN.fa'

        if not os.path.isfile(file_design):
            sys.exit("Design file is needed\n")
        if not os.path.isfile(self.file_out):
            final_out = 'AGE_output.txt'

        PR = PslResult
        #  step 1 load design file
        DESIGNFAS = open(design_fas, 'w')
        design_hash = {}
        with open(file_design, "r") as DESIGN:
            next(DESIGN)
            for line in DESIGN:
                line = line.strip().split("\t")
                id = line[0]
                seq = line[1]
                id = id.replace(' ', '_')
                seq = seq.replace(' ', '')
                line = ">" + id + "\n" + seq + "\n"
                design_hash[id] = seq
                DESIGNFAS.write(line)
        DESIGNFAS.close()

        #######################################################
        #  step 2 - load read files  #############

        # read fastq files
        list_f1 = glob.glob(os.path.join(dat_dir, '*.fq'))
        list_f2 = glob.glob(os.path.join(dat_dir, '*.fastq'))
        list_f1.extend(list_f2)

        ##
        total_read_num = {}
        num_c = 0
        fasta_files = []

        for file in list_f1:
            file_in = file
            dir2, file_short = os.path.split(file)
            fasta_out = file_short + '_INPUT.fa'
            fasta_files.append(fasta_out)

            check_fastq2fasta = PR.fastq2fasta(file_in, fasta_out)
            total_num = check_fastq2fasta
            total_read_num[fasta_out] = total_num

        # read fasta files
        list_f1 = glob.glob(os.path.join(dat_dir, '*.fas'))
        list_f2 = glob.glob(os.path.join(dat_dir, '*.fa'))
        list_f3 = glob.glob(os.path.join(dat_dir, '*.fasta'))
        list_f1.extend(list_f2)
        list_f1.extend(list_f3)

        for file in list_f1:
            fasta_out = file + '_INPUT.fa'
            fasta_files.append(fasta_out)
            total_num = PR.fasta2fasta(file, fasta_out)
            # load number
            total_read_num[fasta_out] = total_num

        ########################################################

        for file_blat_in in fasta_files:
            #  step 3 - blat  ###########################
            blat_out = file_blat_in + '_blat_crispr.psl'
            command_content = ['blat.exe',
                               ' -tileSize=7 -oneOff=1  -maxGap=20 -minIdentity=70 -minScore=20 ',
                               file_blat_in, design_fas, blat_out]
            command_line = " ".join(command_content)
            os.system(command_line)
            print("blat job $file_blat_in is done\n")
            #  step 4 - convert psl -> bed  ##############
            bed_hash_address = PR.psl2bed(blat_out)  # address of one hash
            # step 5 - get sequences number from bed ##############
            get_buf = PR.get_out_buffer
            write_buf = PR.write_buffer_out
            get_ed_note = PR.get_editing_note
            get_ali = PR.get_alignment
            tri_seq = PR.treat_sequence
            tri_in = PR.treat_inter
            reverse = PR.reverseComplement
            get_rev = PR.get_reverse

            read_num_address =  PR.MAIN (file_blat_in, bed_hash_address,self.file_out,get_buf,write_buf,total_read_num,
                                       design_hash,get_ed_note,get_ali,tri_seq,mismatch_cutoff,tri_in,min_cutoff,
                                       indel_report,wt_like_report, get_rev, reverse) # address of one hash

            if remove_files == 1:
                command_line = remove + ' ' + file_blat_in
                os.system(command_line)
                command_line = remove + ' ' + blat_out
                os.system(command_line)
        if remove_files == 1:
            command_line2 = remove + ' ' + design_fas
            os.system(command_line2)
        out += "All Job Done"+"\n"
        self.run_output.setText(out)
        print("All Job Done")


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
# auto


class PslResult(AutoVivification):
    @staticmethod
    def fasta2fasta(input_file, out_file):
        Out_stream = open(out_file, 'w')
        with open(input_file, 'r') as fst:
            all_lines = fst.readlines()
            processed = 0
            for line in "".join(all_lines).split(">")[1:]:
                line = line.split("\n")
                id_out = line[0].split()[0]  # ID actually
                id_out = id_out + '_1'  # add counting number
                seq_out = "".join(line[1:])
                Out_stream.write('>' + id_out + "\n" + seq_out + "\n")
                processed += 1
        return processed

    @staticmethod
    def fastq2fasta(input_file, out_file):
        In_stream = open(input_file, 'r')
        Out_stream = open(out_file, 'w')
        ##
        total_line = 0
        processed = 0
        read_hash = AutoVivification()
        # for loop
        for aline in In_stream:
            aline = aline.rstrip()
            total_line += 1
            if total_line % 4 == 1:
                processed += 1
                id = aline.split()[0][1:]
            elif total_line % 4 == 2:
                seq = aline
                read_hash[seq][id] = 1
        # write Out_stream
        count = 0
        for seq_out in (sorted(read_hash.keys())):
            count += 1
            hits = read_hash[seq_out].keys()
            num = len(hits)
            id_out = "R" + str(count) + '_' + str(num)
            Out_stream.write('>' + id_out + "\n" + seq_out + "\n")
        In_stream.close()
        Out_stream.close()
        return processed

    @staticmethod
    def psl2bed(psl_file):
        psl_array = []
        with open(psl_file, "r") as FILEINBLAT:
            for line in FILEINBLAT.readlines()[5:]:
                lineIn = line
                lineArr = lineIn.strip().split("\t")
                q_size = lineArr[10]
                m_size = lineArr[0]  # match_size
                if int(m_size) / int(q_size) > 0.20:
                    psl_array.append(lineArr)
        psl_line_num = 0
        psl_hash = AutoVivification()
        for lineArr in psl_array:
            psl_line_num += 1
            read_target = lineArr[13]
            query = lineArr[9]
            lineArr = lineArr[8:]
            if not (read_target in psl_hash.keys() and query in psl_hash[read_target].keys()):
                psl_hash[read_target][query] = lineArr
        print(str(psl_line_num) + "total psl lines\n")
        return psl_hash

    @staticmethod
    def reverseComplement(sequence):
        complement = string.maketrans('ATCGN', 'TAGCN')
        return sequence.upper().translate(complement)[::-1]

    @staticmethod
    def treat_inter(q, t):
        q_len = len(q)
        t_len = len(t)
        dis = abs(q_len - t_len)
        oooo = []
        for i in range(dis):
            oooo.append("-")
        if (q_len - t_len) > 0:
            t = t + "".join(oooo)
        else:
            q = q + "".join(oooo)
        return (q, t)

    @staticmethod
    def treat_sequence(q, t, position):
        q_len = len(q)
        t_len = len(t)
        dis = abs(q_len - t_len)

        if dis > 0:
            my_oooo = []
            for i in range(dis):
                my_oooo.append("-")
            # print(my_oooo)
            if q_len > t_len:
                if position == "before":
                    t = "".join(my_oooo) + t
                else:
                    t = t + "".join(my_oooo)
            else:
                if position == "before":
                    q = "".join(my_oooo) + "".join(q)

                else:
                    q = "".join(q) + "".join(my_oooo)
        return (q, t)

    @staticmethod
    def get_reverse(parameter_in, reverseComplement):
        q_string, t_string, indel_string = parameter_in
        q_string = reverseComplement(q_string)
        t_string = reverseComplement(t_string)
        string_test = q_string
        string_test.replace('-', '')
        len_temp = len(string_test)
        indel_in = indel_string.split()
        indel_out = []
        for indel_each in indel_in:
            re_match = re.search(r'^(\d+)(.)(\d+)', indel_each)
            if re_match:
                ooo = str(len_temp - int(re_match.group(1)) - int(re_match.group(3))) + re_match.group(2) + re_match.group(3)
                if re_match.group(2) == 'I':
                    ooo = str(len_temp - int(re_match.group(1))) + re_match.group(2) + re_match.group(3)
                indel_out.insert(0, ooo)
        indel_string = ' '.join(indel_out)
        out_reverse = [q_string, t_string, indel_string]
        return out_reverse

    ##get_reverse

    @staticmethod
    def get_alignment(blocks, blockLengths, qStarts, tStarts, query_seq, target_seq, treat_sequence, treat_inter):

        blockSizes = blockLengths.split(",")
        qArray = qStarts.split(",")
        tArray = tStarts.split(",")
        if blockSizes[-1] == '':
            blockSizes = blockSizes[0:-1]
        if qArray[-1] == '':
            qArray = qArray[0:-1]
        if tArray[-1] == '':
            tArray = tArray[0:-1]

            # target: part of read
        # query : original sequence in design file

        # before blocks
        before_q = query_seq[0:int(qArray[0])]
        before_t = target_seq[0:int(tArray[0])]

        (before_q, before_t) = treat_sequence(before_q, before_t, "before")

        # after blocks
        after_q = query_seq[int(qArray[-1]) + int(blockSizes[-1]):]
        after_t = target_seq[int(tArray[-1]) + int(blockSizes[-1]):]
        (after_q, after_t) = treat_sequence(after_q, after_t, "after")

        # blocks
        med_q = []
        med_t = []
        indel_out = []
        out = ''
        # first block
        med_q_seq = query_seq[int(qArray[0]):(int(blockSizes[0]) + int(qArray[0]))]
        med_t_seq = target_seq[int(tArray[0]):(int(blockSizes[0]) + int(tArray[0]))]

        med_q.append(med_q_seq)
        med_t.append(med_t_seq)
        if int(blocks) > 1:
            for i in range(int(blocks) - 1):
                start_i = qArray[i] + blockSizes[i]
                end_i = int(start_i) + int(qArray[i + 1]) - (int(qArray[i]) + int(blockSizes[i]))
                inter_q_seq = query_seq[int(start_i):int(end_i)]
                start_i = tArray[i] + blockSizes[i]
                end_i = int(start_i) + int(tArray[i + 1]) - (int(tArray[i]) + int(blockSizes[i]))
                inter_t_seq = target_seq[int(start_i):int(end_i)]

                # count deletion insertion
                if len(inter_t_seq) - len(inter_q_seq) > 0:
                    # to fix the mismatch near deletion   length($inter_q_seq)
                    out = str(int(qArray[i]) + int(blockSizes[i]) + len(inter_q_seq)) + 'I' + str(
                        len(inter_t_seq) - len(inter_q_seq))
                if len(inter_q_seq) - len(inter_t_seq) > 0:
                    out = str(int(qArray[i]) + int(blockSizes[i]) + len(inter_t_seq)) + 'D' + str(
                        len(inter_q_seq) - len(inter_t_seq))
                indel_out.append(out)

                (inter_q_seq, inter_t_seq) = treat_inter(inter_q_seq, inter_t_seq)

                med_q.append(inter_q_seq)
                med_t.append(inter_t_seq)

                #  block   after interval
                med_q_seq = query_seq[int(qArray[i + 1]):(int(qArray[i + 1]) + int(blockSizes[i + 1]))]
                med_t_seq = target_seq[int(tArray[i + 1]):(int(tArray[i + 1]) + int(blockSizes[i + 1]))]
                med_q_seq = "".join(med_q_seq)
                med_q.append(med_q_seq)
                med_t.append(med_t_seq)

            # if block
        med_q.append(after_q)
        med_t.append(after_t)

        med_q.insert(0, before_q)
        # print(before_q)
        med_t.insert(0, before_t)

        q_string = ''.join(med_q)
        # print(med_q)
        t_string = ''.join(med_t)

        re_match = re.search(r'(^-+)', q_string)

        if re_match:
            len_temp = len(re_match.group(1))
            q_string = q_string[len_temp:]
            t_string = t_string[len_temp:]

        re_match2 = re.search(r'(-+$)', q_string)
        if re_match2:
            len_temp = len(re_match2.group(1))
            q_string = q_string[0:(-1) * len_temp]

            t_string = t_string[0:(-1) * len_temp]

        out_alignment = [q_string, t_string, ' '.join(indel_out)]
        return out_alignment

    #   get_alignment

    @staticmethod
    def get_editing_note(address_in, query_seq, target_seq, full_read, get_alignment, tr_seq, tri_in, get_reverse,
                          reverseComplement):
        # target: part of read
        # query : original sequence in design file
        (strand, q_name, q_size, q_start, q_end, t_name, t_size, t_start, t_end,
         block_num, block_size, q_start_s, t_start_s) = address_in
        out = []
        if strand == "+":
            out = get_alignment(block_num, block_size, q_start_s, t_start_s, query_seq, full_read, tr_seq, tri_in)
        else:
            out = get_alignment(block_num, block_size, q_start_s, t_start_s, reverseComplement(query_seq), full_read,
                                tr_seq, tri_in)
            out = get_reverse(out, reverseComplement)
        out[0] = "B" + out[0] + "E"
        out[1] = "B" + out[1] + "E"

        if not re.search(r"\w", out[2]):
            out[2] = "no_indel"
        return "\t".join(out)

    ## get_editing_note

    @staticmethod
    def get_out_buffer(fas_file_in, query_with_seq2num_in, seq2alignment_in, total_read_num, min_cutoff, design_hash,
                       indel_report, wt_like_report):

        buffer_out = AutoVivification()
        buffer_num = AutoVivification()
        total_non_redun = 0
        # push output in buffer before printing
        test_case_num = total_read_num[fas_file_in]

        for ref_name in sorted(query_with_seq2num_in.keys()):

            in_hash = query_with_seq2num_in[ref_name]

            # each group
            report_wt_count = 0
            report_indel_count = 0
            other_num = 0
            sub_hit_num = 0
            indel_hit_num = 0
            other_record = []
            # loop for records
            in_hash_sorted = sorted(in_hash.keys())
            for part_of_read in in_hash_sorted:
                hitnum = in_hash[part_of_read]
                hitnum = int(hitnum)

                if hitnum < min_cutoff:
                    continue
                if hitnum / test_case_num < 0.00005:
                    continue
                # get Editing Note for output only
                if part_of_read in seq2alignment_in.keys():
                    (ref_ori_seq, edit_note, array_snp_addr) = seq2alignment_in[part_of_read]
                    (align_ref, align_read, note) = edit_note.split("\t")
                    snp_out = array_snp_addr
                    snp_num = len(snp_out)
                    # deep sequencing SNP
                    if test_case_num > 500 and hitnum <= 2 and snp_num / len(ref_ori_seq) > 0.05:
                        continue
                    ##
                    snp_note = ''
                    if snp_num > 0:
                        snp_note = snp_num + ' SNP(' + ','.join(snp_out)
                    sub_hit_num = sub_hit_num + hitnum
                    total_non_redun += hitnum
                    # output here
                    if ref_ori_seq != design_hash[ref_name]:
                        print("error in read assingment \n")
                    out_line = (
                        fas_file_in, ref_name, ref_ori_seq, part_of_read, hitnum, align_ref, align_read, note, snp_note)
                    # indel
                    if not re.search(r'no.+indel', edit_note):
                        if ((test_case_num > 500 and hitnum <= 3) or (hitnum / test_case_num < 0.001)) and re.search(
                                r'strange', edit_note):
                            continue  # filter strange case in deep sequencing with few reads
                        report_indel_count += 1
                        # total indel
                        indel_hit_num = indel_hit_num + hitnum
                        if report_indel_count <= indel_report:
                            # print REPORT join("\t",@out_line),"\n";
                            # in buffer
                            if not buffer_out[fas_file_in][ref_name]["data"]:
                                bbb = [out_line, ]
                                buffer_out[fas_file_in][ref_name]["data"] = bbb
                            else:
                                buffer_out[fas_file_in][ref_name]["data"].append(out_line)
                        else:
                            other_num = other_num + hitnum
                            other_record = [fas_file_in, ref_name, ref_ori_seq, 'others', hitnum, '-', '-', '-', '-']
                    # no + indel
                    else:  # wt like report
                        report_wt_count += 1
                        if report_wt_count <= wt_like_report:
                            # print REPORT join("\t",@out_line),"\n";
                            # in buffer
                            if not buffer_out[fas_file_in][ref_name]["data"]:
                                bbb = [out_line, ]
                                buffer_out[fas_file_in][ref_name]["data"] = bbb
                            else:
                                buffer_out[fas_file_in][ref_name]["data"].append(out_line)
                        else:
                            other_num = other_num + hitnum
                            other_record = [fas_file_in, ref_name, ref_ori_seq, 'others', hitnum, '-', '-', '-', '-']
                    # wt like report
                # if in key
            # for part
            # other report
            if other_num > 0:
                other_record[-5] = other_num
                # print REPORT  join("\t", @ other_record), "\n";
                if not buffer_out[fas_file_in][ref_name]["data"]:
                    bbb = [other_record, ]
                    buffer_out[fas_file_in][ref_name]["data"] = bbb
                else:
                    buffer_out[fas_file_in][ref_name]["data"].append(other_record)
            # if other num
            # total read number
            total_num = total_read_num[fas_file_in]
            #    0           1          2
            # my @summary_record = ($fas_file_in,$ref_name,'Total Reads: '., 'Total Hits: '.
            # $total_non_redun,'Sub Hits: '.$sub_hit_num, 'Indel Hits: '.$indel_hit_num);
            summary_record = ['Sub Hits: ' + str(sub_hit_num), 'Indel Hits: ' + str(indel_hit_num)]
            buffer_out[fas_file_in][ref_name]['sum'] = summary_record
            buffer_num[fas_file_in]['total'] = total_num
            buffer_num[fas_file_in]['sub'] = total_non_redun
            # print REPORT join("\t",@summary_record),"\n";
        # for ref name
        out_buffer_final = [buffer_out, buffer_num]
        return out_buffer_final

    ##
    @staticmethod
    def write_buffer_out(hash_out, hash_out_num, file):
        print_tracking = 0
        File = open(file, "w")
        for fas_file_in in sorted(hash_out.keys()):
            total_num = hash_out_num[fas_file_in]["total"]
            total_non_redun = hash_out_num[fas_file_in]["sub"]
            for ref_name in sorted(hash_out[fas_file_in].keys()):
                if not hash_out[fas_file_in][ref_name]["data"]:
                    print("No data for $fas_file_in $ref_name \n")
                    continue
                data = hash_out[fas_file_in][ref_name]["data"]
                for da in data:
                    print_tracking += 1
                    if print_tracking == 1:
                        File.write("INPUT\tTarge\tTargetSequence\tReadSequence\tRead#\tAlignedTarget\t"
                                   "AlignedRead\tIndels\tSNPs\n")
                    Da = []
                    for c in da:
                        c = str(c)
                        Da.append(c)
                    File.write("\t".join(Da) + "\n")
                sum_p2 = hash_out[fas_file_in][ref_name]["sum"]
                sum_p1 = (fas_file_in + "\t" + ref_name + "\t" + 'Total Reads:' + str(
                    total_num) + "\t" + 'Total Hits:'
                          + str(total_non_redun))
                File.write(sum_p1 + "\t" + "\t".join(sum_p2) + "\n")

            # print summary region

            File.write("Summary\t -----------------------------------------------------\n")
            File.write("Sum:INPUT\tTarget\tAlignedTarget\tAlignedRead\tTotal Hits\tSub Hits\t"
                       "Indel or WT Hits\tIndel or WT rate %\tPattern\n")

            for ref_name in sorted(hash_out[fas_file_in].keys()):
                if not hash_out[fas_file_in][ref_name]["data"]:
                    continue
                data = hash_out[fas_file_in][ref_name]["data"]

                wt_pair = '' + "\t" + ''
                indel_pair = '' + "\t" + ''
                indel_out = ''
                i = 0
                for da in data:
                    i += 1
                    In = da
                    if not re.search(r"no _indel", In[7]) and In[7] != "-":
                        if indel_pair == "" + "\t" + "":
                            indel_pair = In[5] + "\t" + In[6]
                            indel_temp = In[7].split()

                            ooo_indel = []
                            for ooo_in in indel_temp:
                                re_match1 = re.search(r"I(\d+)", ooo_in)
                                if re_match1:
                                    ooo_indel.append(("+" + re_match1.group(1)))
                                re_match2 = re.search(r"D(\d+)", ooo_in)
                                if re_match2:
                                    ooo_indel.append(("-" + re_match2.group(1)))
                                if re.search(r"strange", ooo_in):
                                    ooo_indel.append("strange editing")
                            indel_out = "Indel:" + ",".join(ooo_indel)
                    else:
                        # no indel
                        if wt_pair != '' + "\t" + '':
                            wt_pair = In[5] + "\t" + In[6]

                # first record

                # sum
                (sub_string, indel_string) = hash_out[fas_file_in][ref_name]["sum"]

                indel_num = 0
                sub_num = 0

                # print 		$sub_string. "sub summary\n";
                # print 		$indel_string ."sub summary\
                re_match3 = re.search(r":[^\d]+(\d+)", sub_string)
                if re_match3:
                    sub_num = re_match3.group(1)
                re_match4 = re.search(r":[^\d]+(\d+)", indel_string)
                if re_match4:
                    indel_num = re_match4.group(1)
                rate = (int(indel_num) / int(sub_num) * 10000) / 100

                out_sum = (fas_file_in, ref_name, indel_pair, total_non_redun, sub_num, indel_num, rate, indel_out)
                File.write('Sum:' + "\t".join("%s" % id for id in out_sum) + "\n")

                wt_num = int(sub_num) - int(indel_num)
                wt_rate = 100 - rate
                out_sum = (fas_file_in, ref_name, wt_pair, total_non_redun, sub_num, wt_num, wt_rate, "WT_like")
                File.write('Sum:' + "\t".join("%s" % id for id in out_sum) + "\n")
            File.write("\n")

    @staticmethod
    def MAIN(fas_file_in, bed_file_address_in, file, get_out_buffer, write_buffer_out, total_read_num, design_hash,
             get_editing_note, get_ali, tri_seq, mismatch_cutoff, tri_in,min_cutoff, indel_report, wt_like_report,
             get_reverse, reverseComplement):
        #########################################################
        # read convereted reads fasta file
        FAS = open(fas_file_in, "r")
        query_with_seq2num = AutoVivification()
        seq2alignment = AutoVivification()
        seq_part2assignment = AutoVivification()
        bed_hash_in = bed_file_address_in
        for lineIn in FAS:
            lineIn = lineIn.rstrip()
            if lineIn.startswith(">"):
                id = lineIn[1:]
                total_num = id.split('_')[-1]
                lineIn = FAS.readline().rstrip()
                # print ("OOOO" + lineIn +"\n")
                seq = lineIn
                if id in bed_hash_in.keys():
                    # psl_hash{$read_target}{$query}  # format of hash
                    # asign read to target reference
                    read_asignment = []
                    # find the best hit in a loop
                    for read_key in bed_hash_in[id].keys():
                        # get annotation from array
                        # Q: design T: read
                        get_array = bed_hash_in[id][read_key]
                        (strand, q_name, q_size, q_start, q_end,
                         t_name, t_size, t_start, t_end,
                         block_num, block_size, q_start_s, t_start_s) = get_array
                        seq_out = seq[int(t_start):int(t_end)]
                        if strand == '-':
                            seq_out = reverseComplement(seq_out)
                        # get original sequence
                        ref_ori_seq = ''
                        if q_name in design_hash.keys():
                            ref_ori_seq = design_hash[q_name]

                        edit_note = get_editing_note(get_array, ref_ori_seq, seq_out, seq,get_ali,tri_seq,tri_in,
                                                     get_reverse,reverseComplement)

                        (align_ref, align_read, note) = edit_note.split("\t")
                        len_ori = len(ref_ori_seq)
                        align_1 = align_ref
                        align_2 = align_read
                        site_snp_list = []
                        iden_num = 0

                        for idx in range(len(align_1)):
                            a1 = align_1[idx]
                            a2 = align_2[idx]
                            if a1 != '-' and a2 != '-':
                                if a1 != a2:
                                    site_snp_list.append(idx)
                            if a1 == a2:
                                iden_num += 1
                        iden_num -= 2
                        # filter on SNP number

                        snp_num_for_filter = len(site_snp_list)

                        dis_max = 0
                        dis_current = 0
                        continuous_site_temp = []
                        continuous_site = []
                        site_snp_checking = site_snp_list

                        if site_snp_checking == []:
                            start = 0
                        else:
                            start = site_snp_checking[0]
                        site_snp_checking = site_snp_checking[1:]
                        for site_checking in site_snp_checking:
                            dis_in = site_checking - start
                            if dis_in < 4:
                                dis_current += 1
                                continuous_site_temp.append(site_checking)
                                if dis_current > dis_max:
                                    dis_max = dis_current
                                    continuous_site = continuous_site_temp
                            else:
                                dis_current = 0
                                continuous_site_temp = []
                            start = site_checking

                        if snp_num_for_filter / len_ori > 0.5:  # first filtering
                            if dis_max < 5:
                                continue

                        if dis_max >= 5:
                            site_snp = []
                            align_x = list(align_read)

                            for iii in continuous_site:
                                align_x[iii] = 'N'

                            align_read = ''.join(align_x)
                            note = 'strange editing'
                            temp_tulp = [align_ref, align_read, note]
                            edit_note = "\t".join(temp_tulp)
                        out_temp = [iden_num, ref_ori_seq, seq_out,
                                    edit_note, site_snp_list, get_array]
                        read_asignment.append(out_temp)
                    # loop for best hit
                    if len(read_asignment) < 1:
                        continue
                    #  assign the first hit
                    read_asignment = sorted(read_asignment, key=lambda x: x[0], reverse=True)
                    (iden_num, ref_ori_seq, seq_out, edit_note,
                     array_snp_addr, bed_array_addr) = read_asignment[0]
                    snp_filtering = array_snp_addr
                    # print ref_ori_seq

                    if len(snp_filtering) / len(ref_ori_seq) > mismatch_cutoff:
                        continue
                    # after assignment
                    # check whether seq_out is assigned earlier
                    # $seq_out has been asigned to one record
                    # just assign reads to that one althogh it may be get one new assignment
                    if seq_out in seq_part2assignment.keys():
                        first_assignment = seq_part2assignment[seq_out]
                        query_with_seq2num[first_assignment][seq_out] += int(total_num)
                        continue
                    # number for part of reads
                    (strand, q_name, q_size, q_start, q_end,
                     t_name, t_size, t_start, t_end,
                     block_num, block_size, q_start_s, t_start_s) = bed_array_addr
                    # this the first asignment
                    if seq_out not in seq2alignment:
                        seq2alignment[seq_out] = [ref_ori_seq, edit_note, array_snp_addr]
                    # assigned
                    if seq_out not in seq_part2assignment:
                        seq_part2assignment[seq_out] = q_name
                        query_with_seq2num[q_name][seq_out] = int(total_num)

        ##
        FAS.close()
        # summary data for output

        (hash_addr1, hash_addr2) = get_out_buffer(fas_file_in, query_with_seq2num, seq2alignment, total_read_num,
                                                  min_cutoff, design_hash, indel_report, wt_like_report)

        # write output
        write_buffer_out(hash_addr1, hash_addr2,file)


##########################################################
# step 1 load design file
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MainUi()
    window.show()
    sys.exit(app.exec_())

