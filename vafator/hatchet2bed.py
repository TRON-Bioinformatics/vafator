import pandas as pd
import pyBigWig


def run_hatchet2bed(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    cn_columns = sorted(list(filter(lambda c: c.startswith('cn_clone'), df.columns)))
    u_columns = sorted(list(filter(lambda c: c.startswith('u_clone'), df.columns)))

    with pyBigWig.open(output_file, "w") as bw:

        chromosomes = []
        starts = []
        ends = []
        values = []
        for i, row in df.iterrows():
            numerator = []
            total_u = 0
            for cn_column, u_column in zip(cn_columns, u_columns):
                u = float(row[u_column])
                total_u += u
                cn = sum(map(lambda c: float(c), row[cn_column].split('|')))
                numerator.append(u * cn)
            chromosomes.append(row['#CHR'])
            starts.append(int(row['START']))
            ends.append(int(row['END']))
            values.append(sum(numerator)/total_u)

            if len(chromosomes) >= 100:
                bw.addEntries(chromosomes, starts, ends=ends, values=values)
                chromosomes = []
                starts = []
                ends = []
                values = []
                break

        if len(chromosomes) > 0:
            bw.addEntries(chromosomes, starts, ends=ends, values=values)
