from numpy.core.umath import log2

__author__ = 'chernov$sergey'

import numpy as np
import pandas as pd
from pandas.io.excel import ExcelWriter, _XlsxWriter
import os
import re
from pandas.io.pytables import HDFStore
import fnmatch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from StringIO import StringIO


def extract_dirlist(_path_):
    listdir = [_path_ + x for x in os.listdir(_path_)]
    listdir = filter(os.path.isdir, listdir)
    return listdir


def extract_filelist(_path_, _extension_):
    matches = []
    for root, dirnames, filenames in os.walk(_path_):
        for filename in fnmatch.filter(filenames, _extension_):
            matches.append(os.path.join(root, filename))
    return matches


def parse_one_and_save(input_file, output_store_name):
    sheet_name = 'All sites'
    skip_rows = [0]
    store = HDFStore(output_store_name)
    df = pd.ExcelFile(input_file).parse(sheetname=sheet_name,
                                        skiprows=skip_rows)
    name = (input_file.split('/')[1]).split('.')[0]
    print "Parsing ", name
    store[name] = df
    store.close()


def parse_list_and_save(list_of_files, output_store_name):
    sheet_name = 'All sites'
    skip_rows = [0]
    store = HDFStore(output_store_name)
    for _file_ in list_of_files:
        df = pd.ExcelFile(_file_).parse(sheetname=sheet_name,
                                        skiprows=skip_rows)
        name = (_file_.split('/')[2]).split('.')[0]
        print "Parsing ", name
        store[name] = df
    store.close()


def clear_sequence(s):
    s = re.sub(r'[_][0-9]', '', s)
    digits = re.findall(r"[-+]?\d*\.\d+|\d+", s)
    for i in digits:
        if float(i) < 0.75:
            s = re.sub(r'\(%s\)' % i, '', s)
        else:
            s = re.sub(i, '1', s)
    return s


def filter_columns(storename, dataframe, collist):
    store = HDFStore(storename)
    df = store[dataframe]
    df = df[collist]
    store[dataframe] = df
    store.close()


def drop_with_low_probability(storename, df_name, loc_probability_colname, threshold=0.95):
    print 'Filtering by low probability in', df_name
    store = HDFStore(storename)
    df = store[df_name]
    if loc_probability_colname is not None:
        df = df[df[loc_probability_colname] >= threshold]
    store[df_name] = df
    store.close()


def fetch_identity_from_web(seq, database='swissprot'):
    """

    :param seq:
    :param database:
    :return: Returns identical proteins in HUMAN and MOUSE organisms
    """""

    def extract_prot_id(string):
        return string.split('|')[3]

    result = []

    record = SeqRecord(Seq(seq), id="tmp", name="", description="")
    SeqIO.write(record, "tmp.fastaa", "fasta")

    # first get the sequence we want to parse from a FASTA file
    f_record = next(SeqIO.parse('tmp.fastaa', 'fasta'))

    # print('Doing the BLAST and retrieving the results...')
    result_handle = NCBIWWW.qblast('blastp', database, f_record.format('fasta'))

    # save the results for later, in case we want to look at it
    save_file = open('my_output.xml', 'w')
    blast_results = result_handle.read()
    save_file.write(blast_results)
    save_file.close()

    # print('Parsing the results and extracting info...')

    # option 1 -- open the saved file to parse it
    # option 2 -- create a handle from the string and parse it
    string_result_handle = StringIO(blast_results)
    b_record = NCBIXML.read(string_result_handle)

    for alignment in b_record.alignments:
        if ('HUMAN' in alignment.title) or ('MOUSE' in alignment.title):
            for hsp in alignment.hsps:
                if hsp.positives == hsp.identities:
                    result.append(extract_prot_id(alignment.title))

    print seq, result, '\n'
    return result


def fetch_indentity_from_local(seq):
    def extract_prot_id(string):
        s = string.split('|')[2]
        s = s.split(' ')[1]
        return s

    result = []
    record = SeqRecord(Seq(seq), id="tmp", name="", description="")
    SeqIO.write(record, "tmp.fastaa", "fasta")

    NcbiblastpCommandline(query='tmp.fastaa', db='_data_/_db_/HUMAN_DB', outfmt=5, out='blastp_human_output.xml')()
    NcbiblastpCommandline(query='tmp.fastaa', db='_data_/_db_/RODENTS_DB', outfmt=5, out='blastp_rodents_output.xml')()

    result_handle = open("blastp_human_output.xml")
    b_record = NCBIXML.read(result_handle)
    for alignment in b_record.alignments:
        for hsp in alignment.hsps:
            if hsp.positives == hsp.identities:
                result.append(extract_prot_id(alignment.title))

    result_handle = open("blastp_rodents_output.xml")
    b_record = NCBIXML.read(result_handle)
    for alignment in b_record.alignments:
        for hsp in alignment.hsps:
            if hsp.positives == hsp.identities:
                result.append(extract_prot_id(alignment.title))

    return ";".join(result)


def fetch_indentity_from_local_batch(seq_list):
    def extract_prot_id(string):
        s = string.split('|')[2]
        s = s.split(' ')[1]
        return s

    LEN = len(seq_list)
    sequences = [0] * LEN
    result_human = []
    result_mouse = []
    for i in np.arange(LEN):
        sequences[i] = SeqRecord(Seq(seq_list[i]), id="seq-" + str(i), description="")

    output_handle = open("all_seq.fasta", "w")
    SeqIO.write(sequences, output_handle, "fasta")
    output_handle.close()

    NcbiblastpCommandline(query='all_seq.fasta', db='_data_/_db_/HUMAN_DB', outfmt=5,
                          out='all_seq_human.xml', num_threads=2)()

    NcbiblastpCommandline(query='all_seq.fasta', db='_data_/_db_/RODENTS_DB', outfmt=5,
                          out='all_seq_rodents.xml', num_threads=2)()

    for record in NCBIXML.parse(open("all_seq_human.xml")):
        if record.alignments:
            tmp = []
            for align in record.alignments:
                for hsp in align.hsps:
                    if hsp.positives == hsp.identities:
                        gene = extract_prot_id(align.title)
                        if len(gene) == 0:
                            print 'oops, extracted empty gene name', align.title
                        else:
                            tmp.append(gene)
            if len(tmp) != 0:
                result_human.append(';'.join(tmp))
            else:
                result_human.append('')
        else:
            result_human.append('')

    for record in NCBIXML.parse(open("all_seq_rodents.xml")):
        if record.alignments:
            tmp = []
            for align in record.alignments:
                for hsp in align.hsps:
                    if hsp.positives == hsp.identities:
                        gene = extract_prot_id(align.title)
                        if len(gene) == 0:
                            print 'oops, extracted empty gene name', align.title
                        else:
                            tmp.append(gene)
            if len(tmp) != 0:
                result_mouse.append(';'.join(tmp))
            else:
                result_mouse.append('')
        else:
            result_mouse.append('')

    print "Length is equivalent? ", len(result_human) == len(sequences) == len(result_mouse)
    return result_human, result_mouse


def make_summary(newcols):
    """

    :param newcols: column names in the main summary table
    :return: none
    """
    print "Making summary..."

    # open store end read base dataframe
    store = HDFStore('_data_/ProteinDataStore.h5')
    df1 = store['Mol_Cell_Proteomics_2011_Epub_2011_September1Supp2']

    # clean sequences
    LEN = len(df1)
    positions = [0] * LEN
    real_glygly = [0] * LEN
    clean_glygly = [0] * LEN
    for i in np.arange(LEN):
        positions[i] = df1['Position'].values[i]
        real_glygly[i] = clear_sequence(df1['GlyGly (K) Probabilities'].values[i])
        clean_glygly[i] = re.sub(r'[^A-Z]', '', real_glygly[i])

    # align with SwissProt Human and Rodents using blastp
    blastpID_HUMAN, blastpID_RODENTS = fetch_indentity_from_local_batch(clean_glygly)

    del df1
    print "Length test", len(positions) == len(real_glygly) == len(clean_glygly) == len(blastpID_HUMAN) == len(
        blastpID_RODENTS)

    # convert to pandas series
    clean_glygly = pd.Series(clean_glygly)
    blastpID_HUMAN = pd.Series(blastpID_HUMAN)
    blastpID_RODENTS = pd.Series(blastpID_RODENTS)

    # Create empty dataframe
    data_summary = pd.DataFrame(columns=newcols)

    # Combine everything required in dataframe
    data_summary['Position'] = positions
    data_summary['GlyGly (K) Probabilities'] = real_glygly
    data_summary['GlyGly Probabilities'] = clean_glygly
    data_summary['SP_ID_BLASTP_HUMAN'] = blastpID_HUMAN
    data_summary['SP_ID_BLASTP_RODENTS'] = blastpID_RODENTS

    # Save to HDF store
    store['DataBases_Summary'] = data_summary
    store.close()


def append_main_summary(newcols, storename_to_append, new_store_position_colname, new_store_glyglyseq_colname):
    """
    This functions only appends data into summary dataframe.
    Data - comliant columns to main dataframe but with
    sequences which are not present  in main dataframe

    :param newcols: column names in the main summary
    :param storename_to_append: store name to open and analyze
    :param new_store_position_colname: Position values column name in the new store
    :param new_store_glyglyseq_colname: Sequences values column name in the new store
    :return: appends main summary table with data, saves changes in the same dataframe in HDF store
    """
    # Open summary from hdf
    store = HDFStore('_data_/ProteinDataStore.h5')
    data_summary = store['DataBases_Summary']
    what_to_append = store[storename_to_append]

    # Get GlyGly values
    all_glygly_in_newstore = what_to_append[new_store_glyglyseq_colname].values
    all_glygly_in_summary = data_summary['GlyGly (K) Probabilities'].values

    # Find new sequences to add to the summary
    newcomer_seqs = []
    for x in all_glygly_in_newstore:
        clear_seq = clear_sequence(x)
        if clear_seq not in all_glygly_in_summary:
            newcomer_seqs.append(clear_seq)

    print len(newcomer_seqs), ' new sequences were found in ', storename_to_append

    # Clean them as well
    clean_newcomer_seqs = map(lambda x: re.sub(r'[^A-Z]', '', x), newcomer_seqs)

    # Find comliant positions
    subset_index = what_to_append[new_store_glyglyseq_colname].isin(
        newcomer_seqs)  # fetch subset where newcomers presents
    positions = what_to_append[subset_index][new_store_position_colname].values

    # BlastP query results
    blastpID_HUMAN, blastpID_RODENTS = fetch_indentity_from_local_batch(clean_newcomer_seqs)

    # convert to pandas series
    positions = pd.Series(positions)
    newcomer_seqs = pd.Series(newcomer_seqs)
    clean_newcomer_seqs = pd.Series(clean_newcomer_seqs)
    blastpID_HUMAN = pd.Series(blastpID_HUMAN)
    blastpID_RODENTS = pd.Series(blastpID_RODENTS)

    # Create empty dataframe to be appended
    data_summary_appendix = pd.DataFrame(columns=newcols)

    # Combine everything required in dataframe
    data_summary_appendix['Position'] = positions
    data_summary_appendix['GlyGly (K) Probabilities'] = newcomer_seqs
    data_summary_appendix['GlyGly Probabilities'] = clean_newcomer_seqs
    data_summary_appendix['SP_ID_BLASTP_HUMAN'] = blastpID_HUMAN
    data_summary_appendix['SP_ID_BLASTP_RODENTS'] = blastpID_RODENTS

    # Append main DataBases_Summary
    data_summary = data_summary.append(data_summary_appendix)

    # Save to HDF store
    store['DataBases_Summary'] = data_summary
    store.close()


def analyze_existence(storename_to_append, gly_gly_seq_colname):
    print "Analyzing occurence in ", storename_to_append

    store = HDFStore('_data_/ProteinDataStore.h5')
    data_summary = store['DataBases_Summary']
    tmp_store_sequences = store[storename_to_append][gly_gly_seq_colname].values
    tmp_store_sequences = map(clear_sequence, tmp_store_sequences)

    # Make binary vector which represents existence
    # of the sequence in storename_to_append dataset
    existense_index = data_summary['GlyGly (K) Probabilities'].isin(tmp_store_sequences).values
    existense_index = np.asarray(existense_index, dtype=int)

    # Create new column in summary table
    data_summary[storename_to_append] = existense_index
    print np.sum(data_summary[storename_to_append])

    # Save to HDF store
    store['DataBases_Summary'] = data_summary
    store.close()


def quantitative_analysis(df_name, df_seq_col, df_quant_col, func=lambda x: x):
    print "Quantitative analysis of ", df_name

    store = HDFStore('_data_/ProteinDataStore.h5')
    summary = store['DataBases_Summary']
    df = store[df_name]
    df = df[[df_seq_col, df_quant_col]]
    renamed_col = '_'.join(df_quant_col.split(' '))
    print "Filling column ", renamed_col
    summary[renamed_col] = ['.'] * len(summary)
    print "Current summary shape: ", summary.shape

    seq_list = map(lambda x: re.sub(r'[^A-Z]', '', x), df[df_seq_col].values)
    for i in zip(seq_list, df[df_quant_col].values):
        query = np.where(summary['GlyGly Probabilities'] == i[0])[0]
        if len(query) != 0:
            index = query[0]
        else:
            print "Omitted data: ", i
            continue

        if not np.isnan(i[1]):
            try:
                tmp = func(i[1])
                summary.loc[index, renamed_col] = tmp
            except Exception as e:
                print i
                print e.message
        else:
            summary.loc[index, renamed_col] = '.'

    store['DataBases_Summary'] = summary
    store.close()


def dump_summary_to_excel(output_filename):
    # Save to XLSX
    store = HDFStore('_data_/ProteinDataStore.h5')
    data_summary = store['DataBases_Summary']
    writer = ExcelWriter(output_filename + '.xlsx', engine='xlsxwriter')
    data_summary.to_excel(writer, 'DataBases_Summary', index=True)
    writer.save()


def colorful_dump_summary_to_excel(output_filename, range_label='L1:U36229'):
    # < -2 dark green
    # -2 to -1 light green
    # -1 to  1 yellow
    # 1 to 2 Orange
    # > 2 red
    store = HDFStore('_data_/ProteinDataStore.h5')
    data_summary = store['DataBases_Summary']
    writer = ExcelWriter(output_filename + '.xlsx', engine='xlsxwriter')
    data_summary.to_excel(writer, 'DataBases_Summary', index=True)

    workbook = writer.book
    worksheet = writer.sheets['DataBases_Summary']

    # using pallete http://www.colourlovers.com/palette/3687876/
    blue = workbook.add_format({'bg_color': '#69D2E7', 'font_color': '#000000'})
    coral = workbook.add_format({'bg_color': '#A7DBD8', 'font_color': '#000000'})
    yellow = workbook.add_format({'bg_color': '#EAE319', 'font_color': '#000000'})
    orange = workbook.add_format({'bg_color': '#FA6900', 'font_color': '#000000'})
    red = workbook.add_format({'bg_color': '#E2434B', 'font_color': '#000000'})
    # empty = workbook.add_format({'bg_color': '#FFFFFF', 'font_color': '#000000'})
    #
    # worksheet.conditional_format(range_label, {'type': 'text',
    #                                            'criteria': 'begins with',
    #                                            'value': '.',
    #                                            'format': empty})

    worksheet.conditional_format(range_label, {'type': 'cell',
                                               'criteria': '<',
                                               'value': -2,
                                               'format': blue})

    worksheet.conditional_format(range_label, {'type': 'cell',
                                               'criteria': 'between',
                                               'minimum': -2,
                                               'maximum': -1,
                                               'format': coral})

    worksheet.conditional_format(range_label, {'type': 'cell',
                                               'criteria': 'between',
                                               'minimum': -1,
                                               'maximum': 1,
                                               'format': yellow})

    worksheet.conditional_format(range_label, {'type': 'cell',
                                               'criteria': 'between',
                                               'minimum': 1,
                                               'maximum': 2,
                                               'format': orange})

    worksheet.conditional_format(range_label, {'type': 'cell',
                                               'criteria': '>',
                                               'value': 2,
                                               'format': red})
    writer.save()
    store.close()


def reindex_summary():
    store = HDFStore('_data_/ProteinDataStore.h5')
    data_summary = store['DataBases_Summary']
    range_index = [x for x in np.arange(len(data_summary))]
    print "Reindexing..."
    data_summary = data_summary.set_index([range_index])
    store['DataBases_Summary'] = data_summary
    store.close()


def del_col_from_summary(colname):
    store = HDFStore('_data_/ProteinDataStore.h5')
    data_summary = store['DataBases_Summary']
    data_summary = data_summary.drop(colname, 1)
    store['DataBases_Summary'] = data_summary
    store.close()


def main():
    raw_files_path = '_data_/_xls_/'

    storename = '_data_/ProteinDataStore.h5'

    summary_name = '_data_/DataBases_Summary'

    newcols = ['Position',
               'GlyGly (K) Probabilities',
               'GlyGly Probabilities',
               'SP_ID_BLASTP_HUMAN',
               'SP_ID_BLASTP_RODENTS']

    dataset_1 = {'name': 'Mol_Cell_Proteomics_2011_Epub_2011_September1Supp2',
                 'prob_col': 'Localization Probability',
                 'pos_col': 'Position',
                 'seq_col': 'GlyGly (K) Probabilities'}

    dataset_2 = {'name': 'Mol_Cell_Proteomics_2011_Epub_2011_September1Supp3',
                 'prob_col': 'Localization Probability',
                 'pos_col': 'Position',
                 'seq_col': 'GlyGly (K) Probabilities',
                 'quant_col_list': ['Ratio H/L Experiment 1',
                                    'Ratio H/L Experiment 2']}

    dataset_3 = {'name': 'Mol_Cell_Proteomics_1578_Supp3',
                 'prob_col': 'Localization probability',
                 'pos_col': 'Position',
                 'seq_col': 'GlyGly Probabilities'}

    dataset_4 = {'name': 'Mol_Cell_Proteomics_148_Supp2_dataset_4',
                 'prob_col': None,
                 'pos_col': 'Site',
                 'seq_col': 'Peptide Identifier',
                 'quant_col_list': ['Rep1 SILAC Log2 Ratio 5uM MG-132/0.5% DMSO (H/L)',
                                    'Rep2 SILAC Log2 Ratio 5uM MG-132/0.5% DMSO (L/M)',
                                    'Rep1 SILAC Log2 Ratio 5uM PR-619/0.5% DMSO (M/L)',
                                    'Rep2 SILAC Log2 Ratio 5uM PR-619/0.5% DMSO (H/M)',
                                    'Rep1 SILAC Log2 Ratio 17uM PR-619/0.5% DMSO (M/L)',
                                    'Rep2 SILAC Log2 Ratio 17uM PR-619/0.5% DMSO (H/M)',
                                    'Rep1 SILAC Log2 Ratio 17uM PR-619 & 5uM MG-132/0.5% DMSO (H/L)',
                                    'Rep2 SILAC Log2 Ratio 17uM PR-619 & 5uM MG-132/0.5% DMSO (L/M)']}

    dataset_5 = {'name': 'Science_834_SupplementaryTableS1',
                 'prob_col': 'Localization Prob',
                 'pos_col': 'Position',
                 'seq_col': 'Acetyl (K) Probabilities'}

    df_list = [dataset_1, dataset_2, dataset_3, dataset_4, dataset_5]

    ###
    # = = = = = = = = PROCESSING PIPELINE = = = = = = = =
    ###
    # excel_filenames = extract_filelist(raw_files_path, '*.xlsx')
    # parse_list_and_save(excel_filenames, storename)
    #
    # for store in df_list:
    # drop_with_low_probability(storename, store['name'], store['prob_col'])
    #
    # make_summary(newcols)
    #
    # for store in df_list[1:]:
    #     append_main_summary(newcols, store['name'], store['pos_col'], store['seq_col'])
    # reindex_summary()

    # for store in df_list:
    #     analyze_existence(store['name'], store['seq_col'])
    #
    # for i in dataset_2['quant_col_list']:
    #     quantitative_analysis(dataset_2['name'], dataset_2['seq_col'], i, func=lambda x: log2(x))
    #
    # for i in dataset_4['quant_col_list']:
    #     quantitative_analysis(dataset_4['name'], dataset_4['seq_col'], i)

    colorful_dump_summary_to_excel(summary_name)
    # dump_summary_to_excel(summary_name)


if __name__ == '__main__':
    main()