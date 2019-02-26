from datetime import datetime
import os
import hashlib
from copy import deepcopy

from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas
import jinja2
import weasyprint

import flametree
from pdf_reports import (write_report, pug_to_html, dataframe_to_html,
                         style_table_rows, add_css_class)
import pdf_reports

from .version import __version__

THIS_PATH = os.path.dirname(os.path.realpath(__file__))
ASSETS_PATH = os.path.join(THIS_PATH, "reports_assets")
DOMESTICATION_REPORT_TEMPLATE = os.path.join(ASSETS_PATH,
                                             "domestication_report.pug")
STYLESHEET = os.path.join(ASSETS_PATH, 'report_style.css')


def genedom_pug_to_html(template, **context):
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    defaults = {
        'genedom_sidebar_text': "Generated on %s by Genedom version %s" %
                                (now, __version__),
        'genedom_logo_url': os.path.join(ASSETS_PATH, 'imgs', 'logo.png'),
    }
    for k in defaults:
        if k not in context:
            context[k] = defaults[k]
    return pug_to_html(template, **context)


def write_pdf_domestication_report(target, domestication_infos, domesticators):
    summary_table = dataframe_to_html(
        domestication_infos,
        extra_classes=("definition",)
    )
    def tr_modifier(tr):
        tds = list(tr.find_all("td"))
        if len(tds) == 0:
            return
        if len(tds) == 6:
            name, order_id, dom, final, added, edited = tds
        else:
            name, order_id, dom, final, added, edited, barcode = tds
        if "Failed" in final.text:
            add_css_class(tr, "negative")
        if edited.text != "0":
            add_css_class(edited, "warning")
    summary_table = style_table_rows(summary_table, tr_modifier)
    html = genedom_pug_to_html(
        DOMESTICATION_REPORT_TEMPLATE,
        summary_table=summary_table,
        domesticators=domesticators
    )
    write_report(html, target, extra_stylesheets=(STYLESHEET,))
