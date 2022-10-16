import itertools
import os

import fire
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams.update({'font.size': 22})


def generate_dockopt_job_report(dockopt_job_dir_path=".", pdf_path="dockopt_job_report.pdf", opt_results_csv_file_name="dockopt_job_results.csv", enrichment_metric="enrichment_score"):
    opt_results_csv_file_path = os.path.join(dockopt_job_dir_path, opt_results_csv_file_name)
    df = pd.read_csv(opt_results_csv_file_path)

    with PdfPages(pdf_path) as f:
        #
        fig = plt.figure(figsize=(11.0, 8.5))

        image_file_path = os.path.join(dockopt_job_dir_path, "best_retrodock_job", "roc.png")
        image = mpimg.imread(image_file_path)
        plt.axis('off')
        plt.suptitle("linear-log ROC plot of best job")
        plt.imshow(image)

        f.savefig(fig)
        plt.close(fig)

        #
        multivalued_config_param_columns = [column for column in df.columns if "." in column and df[column].nunique() > 1]
        for column in multivalued_config_param_columns:
            fig = plt.figure(figsize=(11.0, 8.5))

            if len(df) == len(df.drop_duplicates(subset=[column])):
                df.plot.bar(x=column, y=enrichment_metric, rot=0)
            else:
                sns.boxplot(data=df, x=column, y=enrichment_metric, showfliers=False, boxprops={'facecolor': 'None'})
                sns.stripplot(data=df, x=column, y=enrichment_metric, zorder=0.5)

            f.savefig(fig)
            plt.close(fig)

        #
        df = df.sort_values(by=enrichment_metric, ascending=False, ignore_index=True)
        for column_1, column_2 in itertools.combinations(multivalued_config_param_columns, 2):
            #
            fig, ax = plt.subplots()
            fig.set_size_inches(11.0, 8.5)

            #
            df_no_duplicates = df.drop_duplicates(subset=[column_1, column_2], keep="first", ignore_index=True)

            #
            df_pivot = pd.pivot_table(df_no_duplicates, values=enrichment_metric, index=[column_1], columns=[column_2])
            df_pivot = df_pivot.sort_index(axis=0, ascending=False)
            sns.heatmap(df_pivot, ax=ax, annot=True, square=True, fmt='.2f', center=0, cmap='icefire', robust=True,
                        cbar_kws={'label': enrichment_metric})

            f.savefig(fig)
            plt.close(fig)


if __name__ == '__main__':
    fire.Fire(generate_dockopt_job_report)
