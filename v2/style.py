import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import re
CAPTIONS={}
def setup():
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,
        'axes.titleweight':'bold','axes.labelsize':8,'xtick.labelsize':7.5,
        'ytick.labelsize':7.5,'legend.fontsize':7.2,'axes.edgecolor':C['ink'],
        'text.color':C['ink'],'axes.labelcolor':C['ink'],'xtick.color':C['ink'],
        'ytick.color':C['ink'],'axes.spines.top':False,'axes.spines.right':False,
        'axes.linewidth':.65,'xtick.major.width':.65,'ytick.major.width':.65,
        'axes.axisbelow':True,'grid.linewidth':.45,'grid.alpha':.35,
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
        'svg.hashsalt':'opv-revision-final','savefig.dpi':600,
        'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic',
        'mathtext.bf':'Arial:bold','mathtext.bfit':'Arial:bold:italic',
        'mathtext.sf':'Arial','mathtext.fallback':'stix'})

def panel(ax, letter, title, pad=10):
    ax.set_title(f'{letter}  {title}',loc='left',pad=pad)

def note(fig, s, x=.08, y=.022, width=120):
    fig.text(x,y,'\n'.join(textwrap.wrap(s,width)),fontsize=7,va='bottom',linespacing=1.35)

def short_label(s):return bool(re.fullmatch(r'(?:P\d+[a-z]?|[AD]\d|SM\d+|\d+)',str(s),re.I))

def save(fig,name,caption):
    for ext in ['svg','pdf','png']:
        fig.savefig(FIG/(name+'.'+ext),facecolor='white',bbox_inches=None,dpi=600)
    CAPTIONS[name]=caption
    plt.close(fig)

import pandas as pd
import helpers as base
WIDTH=180/25.4
def coverage(frame,specs,table):
    rows=[]
    for fields,label,kind in specs:
        mask=frame[list(fields)].apply(base.number).notna().all(axis=1) if kind=='number' else frame[list(fields)].apply(base.present).all(axis=1)
        n=int(mask.sum());rows.append({'table':table,'field_group':label,'fields':';'.join(fields),'nonmissing_records':n,'denominator':len(frame),'percent':100*n/len(frame)})
    return pd.DataFrame(rows)

def coverage_bars(ax,cov,color):
    ax.barh(cov.field_group,cov.percent,color=color,height=.6);ax.invert_yaxis()
    ax.set_xlim(0,132);ax.set_xticks([0,25,50,75,100]);ax.set_xlabel('Records with information (%)');ax.grid(axis='x')
    for i,r in cov.reset_index(drop=True).iterrows():ax.text(103,i,f'{r.nonmissing_records:,}',va='center',fontsize=7.3)
