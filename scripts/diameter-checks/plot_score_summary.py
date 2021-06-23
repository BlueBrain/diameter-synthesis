"""Plotting functions."""

# Copyright (C) 2021  Blue Brain Project, EPFL
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import division

import json
import os
import collections

import matplotlib
matplotlib.use('Agg')

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

CLIP_LIMIT = 10

def plot_scores_summary(plot_dir):

    csv_path = 'output/scores.csv'
    pdf_path = os.path.join(plot_dir, 'scores_test_super_mtype.pdf')
    print('Plotting scores to %s' % pdf_path)

    scores = pd.read_csv(csv_path)

    # Remove original optimization protocol and only keep MM protocols
    # in this case these two are duplicates
    scores = scores.loc[scores['is_original'] == 0]
    exemplar_morphnames = scores.loc[(scores['is_exemplar'] == 1) & (
        scores['is_repaired'] == 1)]['morph_name'].values
    emodel_morphnames = scores.loc[(scores['is_exemplar'] == 1) & (
        scores['is_repaired'] == 1)]['emodel'].values

    n_morph = 100
    with PdfPages(pdf_path) as pdf:
        for exemplar_morphname, exemplar_emodel in zip(
                exemplar_morphnames[:n_morph], emodel_morphnames[:n_morph]):

            this_scores = scores.loc[scores['emodel'].str.contains(
                exemplar_emodel)]

            morph_rows = list(this_scores.iterrows())

            # For every evaluated me-combo, plot the features using a bar plot
            stds = []
            for index, (_, row) in enumerate(morph_rows):

                mecombo_scores = collections.OrderedDict(
                        json.loads(row['scores']))
                mecombo_scores = collections.OrderedDict(
                    sorted(mecombo_scores.items()))

                if row['is_exemplar'] == 1:
                    if row['is_repaired'] == 1:
                        stds_repaired = list(mecombo_scores.values())
                        repaired_text = ' repaired'
                    else:
                        stds_unrepaired = list(mecombo_scores.values())
                else:
                    stds.append(list(mecombo_scores.values()))


            plt.figure(figsize=(10, 10))
            plt.title('%s %s' % (exemplar_morphname, exemplar_emodel))

            stds = np.array(stds)
            stds[stds>CLIP_LIMIT] = CLIP_LIMIT

            stds_repaired = np.array(stds_repaired)
            stds_repaired[stds_repaired>CLIP_LIMIT] = CLIP_LIMIT

            stds_unrepaired = np.array(stds_unrepaired)
            stds_unrepaired[stds_unrepaired>CLIP_LIMIT] = CLIP_LIMIT

            plt.boxplot(stds, vert = False)
            for i, std_rep in enumerate(stds_repaired):
                if std_rep - stds_unrepaired[i]>0:
                    plt.barh(i+1, width = std_rep - stds_unrepaired[i], left = stds_unrepaired[i], color='0.5')
                else:
                    plt.barh(i+1, width = stds_unrepaired[i]- std_rep, left =  std_rep, color='C0')
                width = 1
                plt.axhline(i+0.5,lw=0.5, c='k')

            plt.xlabel('score (#std), clipped at %g' % CLIP_LIMIT)
            #plt.gca().set_yticks(
            #            [t + (len(morph_rows) / 2) * width for t in y])
            plt.gca().set_yticklabels(mecombo_scores.keys(),
                                              rotation=0,
                                              fontsize=6)


            #plt.legend()
            plt.tight_layout()
            pdf.savefig()
            plt.close()

if __name__ == '__main__':
    plot_scores_summary('.')
