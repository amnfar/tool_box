#Amna Farooq 2020
import sys, os
import argparse
import pandas as pd
import numpy as np
import plotly.plotly as py
import plotly
from plotly import tools
import plotly.graph_objs as go
import plotly.io as pio
import time

start_time = time.time()


def read_bedMethyl(in_path):
    bed_methyl = pd.read_csv(in_path, sep='\t', usecols=[0, 1, 2, 3, 4], header=None, comment='#',
                             names=['chr', 'chr_start', 'chr_end', 'meth_perc', 't_reads'],
                             dtype={0: 'str', 1: 'int', 2: 'int', 3: 'float', 4: 'int'})
    bed_methyl['chr'] = bed_methyl.chr.str.replace('chr', '', case=None)

    print "Total elements in input file: ", bed_methyl.shape[0]
    return bed_methyl


def chr_style(df, org, chr, num, d_type):
    pd.set_option('mode.chained_assignment', None) # to suppress copy warning.
    if chr == 'no':
        df['chr'] = df.chr.str.replace('chr', '', case=None)
    if chr == 'yes':
        df['chr'] = 'chr' + df.chr.astype(str)
    if org == "mou":
        if num == 'no':
            df['chr'] = df.chr.str.replace('20', 'X')
            df['chr'] = df.chr.str.replace('21', 'Y')
            df['chr'] = df.chr.str.replace('22', 'M')
        if num == 'yes':
            df['chr'] = df.chr.str.replace('X', '20', case=None)
            df['chr'] = df.chr.str.replace('Y', '21', case=None)
            df['chr'] = df.chr.str.replace('MT', '22', case=None)
            df['chr'] = df.chr.str.replace('M', '22', case=None)

    else:
        if num == 'no':
            df['chr'] = df.chr.str.replace('23', 'X')
            df['chr'] = df.chr.str.replace('24', 'Y')
            df['chr'] = df.chr.str.replace('25', 'M')
        if num == 'yes':
            df['chr'] = df.chr.str.replace('X', '23', case=None)
            df['chr'] = df.chr.str.replace('Y', '24', case=None)
            df['chr'] = df.chr.str.replace('MT', '25', case=None)
            df['chr'] = df.chr.str.replace('M', '25', case=None)

    if d_type == 'str':
        df['chr'] = df.chr.astype('str')
    if d_type == 'num':
        df['chr'] = df.chr.astype('int')
    return df


def get_chr_data(df, chr):
    df_split_chr = df[df.chr == int(chr)]
    return df_split_chr

def clean_sort(df, org, e, out_file):
    df = df[df.chr.str.contains('\.') == False]
    df_num = chr_style(df, org, chr = 'none', num= 'yes',d_type= 'num')
    df_num = df_num.sort_values(by=['chr', 'chr_start'])
    print "Data points left after cleaning = ", df_num.shape[0]
    if e == 'y':
        print "Exporting filtered, sorted version of input."
        print "Number of elements exported = ", df_num.shape[0]
        np.savetxt(out_file, df_num.values, fmt='%s', delimiter='\t')
    return df_num

def avg_read_coverage(df_num, out_ARC):
    
    print "count of meth points",df_num.shape[0]
    total_reads = df_num.t_reads.sum()
    reads_5 = df_num[df_num.t_reads <= 5].count()
    reads_0 = df_num[df_num.t_reads == 0].count()
    
    """reads_10= df_num[df_num.t_reads > 10].count()
        reads_15= df_num[df_num.t_reads > 15].count()
        reads_20 = df_num[df_num.t_reads > 20].count() """
    print "Total sum Reads = " , total_reads
    print "Number of Reads less equal than 5 = ", reads_5[1]
    print "Number of reads equal to 0 = ", reads_0[1]

    step = 5
    intervals = []
    df_ARC = df_num.groupby(pd.cut(df_num["t_reads"], np.arange(0, 21,step))).count()
    for i in np.arange(1,21,step):
        intervals.append(str(i)+"-"+str(i+step-1))
    
    df_ARC_plot = df_ARC ['t_reads']

    df_ARC_plot= df_ARC_plot.multiply(100)
    df_ARC_plot = df_ARC_plot.divide(df_num.shape[0])
    plotly.tools.set_credentials_file(username='DemoAccount', api_key='lr1c37zw81')

    fig = go.Figure(
                data=[go.Bar(y=df_ARC_plot)] ,
                layout = go.Layout(
                                   xaxis = go.layout.XAxis(
                                                           tickmode = 'array',
                                                           tickvals = np.arange(0,len(intervals)),
                                                           ticktext = intervals)))
    fig['layout']['yaxis'].update(title= 'Percentage',range=[0, 100])
    fig['layout']['xaxis'].update(title= 'Number of Reads')
    fig['layout']['title'].update(text= 'Average Read Coverage for '+file_name)

    plotly.offline.plot(fig , filename= out_ARC)
    pio.write_image(fig,  out_ARC , format = 'png')


    return df_num


    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     add_help=True, description="""This script is used for
                                         \n1) Cleaning and sorting input methylation data (removes unidentified chromosmes)
                                         \n2) Extracting Chr wise data from input methylation files.
                                         \n3) Calculates average read count between the range 1-5, 6-10, 11-15, 16-20.
                                         NOTE: Input file must have following order of columns : chr, chr_start, chr_end, meth_perc, t_reads """)
    try:
        required_named = parser.add_argument_group('required arguments')
        
        optional = parser.add_argument_group('optional arguments')
        
        required_named.add_argument("--In_path", help="Enter the path for input file",
                                    required=True, type=str)
        required_named.add_argument("--Out_path", help="Enter the path for output files",
                                    default='out', required=True, type=str)
        required_named.add_argument("--Org", help="Enter hum for human, mou for mouse. Default is hum.",
                                    default='hum', required=True, type=str)
        optional.add_argument("--Chr", help="""Enter the chromosome number to extract data for specific chromosome 
                                    or 'all' to split the input chromosome wise 
                                    or 'none' to treat input file as it is. Default = none""",
                              default = 'none', type=str)
        optional.add_argument("--E", help="Extract sorted version of primary input file [y/n]. Default = n ",
                              default='n', type=str)
        optional.add_argument("--Avg_read",
                              help= " If you want to know average read coverage in the input file, type 'y'.",
                              type=str)
        
        args = parser.parse_args()
    
    except IOError as err:
        
        print >> sys.stderr, "Command line arguments error:", err
        exit(1)
file_name = os.path.basename(args.In_path)


if not os.path.exists(args.Out_path):
    os.makedirs(args.Out_path)

bed_methyl = read_bedMethyl(args.In_path)

if args.Avg_read == 'y':
    out_ARC= args.Out_path + '/ARC_'+file_name + '.png'
    avg_read_coverage(bed_methyl, out_ARC)

if args.E == 'y':
    S_file_name = 'sorted_' + file_name
    out_sorted = args.Out_path + '/' + file_name
    if not os.path.exists(args.Out_path):
        os.makedirs(args.Out_path)
    clean_sort(bed_methyl, args.E , out_sorted)



if args.Chr == 'all':
    sort_bed_methyl = clean_sort(bed_methyl, args.Org, args.E, "none")
    if args.Org == 'hum':
        last_chr = 26
    else:
        last_chr = 23
    for chr in range(1,last_chr):
        get_chr = get_chr_data( sort_bed_methyl, chr)
        print"\n****** Data for CHR"+str(chr)+" ******\n"
        print "Number of elements for Chr" + str(chr) + " exported = " + str(get_chr.shape[0])
        get_chr = chr_style(get_chr,args.Org, chr='yes', num='none', d_type='none')
        file_name1 = "chr"+ str(chr) + '_' + file_name
        path = args.Out_path +'/'+"chr"+ str(chr)
        outfile = path +'/' + file_name1
        if not os.path.exists(path):
            os.makedirs(path)
        np.savetxt(outfile, get_chr.values, fmt='%s', delimiter='\t')


elif args.Chr == 'none':
    sort_bed_methyl = clean_sort(bed_methyl, args.Org, args.E, "none")
    get_chr = chr_style(sort_bed_methyl, args.Org, chr='yes', num='none', d_type='none')
    file_name1 = "chr"  + '_' + file_name
    path = args.Out_path + '/' + "chr"
    outfile = path + '/' + file_name1
    if not os.path.exists(path):
        os.makedirs(path)
    np.savetxt(outfile, get_chr.values, fmt='%s', delimiter='\t')


else:
    sort_bed_methyl = clean_sort(bed_methyl, args.Org, args.E, "none")
    get_chr = get_chr_data( sort_bed_methyl, args.Chr)
    print"\n****** Data for CHR"+str(args.Chr)+" ******\n"
    print "Number of elements for Chr" + str(chr) + " exported = " + str(get_chr.shape[0])
    get_chr = chr_style(get_chr, args.Org, chr='yes', num='none', d_type='none')
    file_name1 = "chr"+ str(args.Chr) + '_' + file_name
    path = args.Out_path +'/'+"chr"+ str(args.Chr)
    outfile =  path +'/' + file_name1
    if not os.path.exists(path):
        os.makedirs(path)
    np.savetxt(outfile, get_chr.values, fmt='%s', delimiter='\t')

print "Time taken(in seconds): ", (time.time() - start_time)
