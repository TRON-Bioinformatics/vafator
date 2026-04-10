import pandas as pd


def run_hatchet2bed(input_file, output_prefix):
    input_df = pd.read_csv(input_file, sep="\t")
    cn_columns = sorted([c for c in input_df.columns if c.startswith("cn_clone")])
    u_columns = sorted([c for c in input_df.columns if c.startswith("u_clone")])

    for sample in input_df.SAMPLE.unique():
        data = []
        for i, row in input_df[input_df.SAMPLE == sample].iterrows():
            numerator = []
            total_u = 0
            for cn_column, u_column in zip(cn_columns, u_columns):
                u = float(row[u_column])
                total_u += u
                cn = sum(map(float, row[cn_column].split("|")))
                numerator.append(u * cn)
            data.append(
                [row["#CHR"], row["START"], row["END"], sum(numerator) / total_u]
            )

        output_df = pd.DataFrame(
            data=data, columns=["chromosome", "start", "end", "value"]
        )
        output_df.to_csv(
            "{}.{}.bed".format(output_prefix, sample),
            sep="\t",
            header=False,
            index=False,
        )
