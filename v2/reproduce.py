from pathlib import Path
import argparse,json,shutil,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import helpers as base
import style
NAMES = ['fig1_workflow_and_scope','fig2_collection_profile',
         'fig3_core_fields_and_performance','fig4_full_data_and_subsets',
         'figS1_supplementary_collections']
C = {'blue':'#365f79','teal':'#4f8584','gray':'#9aa6ae','orange':'#b17c63',
     'purple':'#4f8584','light':'#edf1f3','ink':'#25333c','red':'#a96d55'}
DENSITY_CMAP = LinearSegmentedColormap.from_list('opv_density',['#bad1dc','#83a9bd','#4e7b96','#254a63'])
def csv(frame,name):
    frame.to_csv(SOURCE/name,index=False,encoding='utf-8-sig',float_format='%.8g')

def ranking(ax,frame,labelcol,letter,title,colors=None):
    q=frame.sort_values(['source_dois',labelcol],ascending=[False,True]).head(8).sort_values('source_dois')
    colors=[C['teal'] if base.family(x)=='Fullerene label match' else C['blue'] for x in q[labelcol]] if colors else C['blue']
    ax.barh(q[labelcol],q.source_dois,height=.65,color=colors)
    ax.set_xlim(0,q.source_dois.max()*1.24);ax.grid(axis='x');ax.set_xlabel('Source DOIs')
    ax.xaxis.set_major_locator(MaxNLocator(4,integer=True));ax.tick_params(axis='y',labelsize=7)
    for i,v in enumerate(q.source_dois):ax.text(v+q.source_dois.max()*.025,i,f'{v:,}',va='center',fontsize=6.9)
    style.panel(ax,letter,title)

def vector_colorbar(fig,im,ax,label,ticks=None):
    c=fig.colorbar(im,ax=ax,fraction=.05,pad=.065,ticks=ticks)
    c.solids.set_rasterized(False);c.solids.set_edgecolor('face');c.set_label(label,fontsize=7)
    c.ax.tick_params(labelsize=6.5);return c

def subset_analysis(b):
    """Apply published selection rules to current values, never cached tier flags."""
    rules=json.loads((RELEASE/'validation/subset_rules.json').read_text(encoding='utf-8'))
    ff=b.ff_num.where(b.ff_num.gt(1.5),b.ff_num*100)
    complete=b[['pce_num','voc_num','jsc_num','ff_num']].notna().all(axis=1)
    calculated=b.voc_num*b.jsc_num*ff/100.0
    relative=(calculated-b.pce_num).abs()/b.pce_num.abs()*100
    eligible=complete&b.pce_num.gt(0)
    required=complete&b[['doi_norm','donor','acceptor']].apply(base.present).all(axis=1)
    ranges=b.pce_num.gt(0)&b.pce_num.le(30)&b.voc_num.between(0,2)&b.jsc_num.between(0,50)&ff.between(0,100)
    within2=eligible&relative.le(2)
    candidates=required&ranges&within2
    strict=pd.Series(False,index=b.index);duplicate_of=pd.Series('',index=b.index)
    seen_reported={};seen_canonical={}
    ordered=b.loc[candidates].assign(sort_id=pd.to_numeric(b.loc[candidates,'id'],errors='coerce'),sort_text=b.loc[candidates,'id'].astype(str)).sort_values(['sort_id','sort_text'],na_position='last',kind='stable')
    for idx in ordered.index:
        r=b.loc[idx]
        # Duplicate keys use exact stripped public cell strings.
        core=tuple(str(r[f]).strip() for f in ['doi_norm','voc','jsc','ff','pce'])
        reported=(str(r.donor).strip(),str(r.acceptor).strip())+core
        canon=tuple(str(r[c]).strip() if str(r[c]).strip().lower() not in base.EMPTY else str(r[f]).strip() for c,f in [('donor_canonical','donor'),('acceptor_canonical','acceptor')])+core
        if reported in seen_reported or canon in seen_canonical:
            duplicate_of.at[idx]=seen_reported.get(reported,seen_canonical.get(canon,''));continue
        strict.at[idx]=True;seen_reported[reported]=r.id;seen_canonical[canon]=r.id
    paired=b[['donor_smiles','acceptor_smiles']].apply(base.present).all(axis=1)
    masks={'full':pd.Series(True,index=b.index),'strict_performance':strict,'strict_molecular':strict&paired}
    membership=b[['id','doi']].copy()
    for k,m in [('complete_numeric_quartet',complete),('required_fields_present',required),('within_physical_ranges',ranges),('within_2pct_at_nominal_100',within2),('both_smiles_present',paired),*masks.items()]:membership[k]=m
    membership['duplicate_of_id_under_strict_rule']=duplicate_of
    csv(membership,'figure4_current_subset_membership.csv')
    inputs=b.loc[eligible,['id','doi','voc_num','jsc_num','ff_num','pce_num']].copy()
    inputs['ff_percent_used']=ff[eligible];inputs['nominal_pin_mw_cm2']=100
    inputs['pce_calculated_nominal_percent']=calculated[eligible];inputs['relative_difference_percent']=relative[eligible]
    for key,mask in masks.items():inputs[key]=mask[eligible]
    csv(inputs,'figure4_nominal_one_sun_inputs.csv')
    summaries=[]
    for key,m in masks.items():
        summaries.append({'subset':key,'records':int(m.sum()),'source_dois':int(b.loc[m,'doi_norm'].nunique()),'retained_percent':100*m.mean(),'complete_quartets':int((m&complete).sum()),'paired_smiles':int((m&paired).sum()),'median_reported_pce':b.loc[m,'pce_num'].median()})
    summary=pd.DataFrame(summaries);csv(summary,'figure4_subset_summary.csv')
    threshold=pd.DataFrame([{'maximum_relative_difference_percent':t,'records':int((eligible&relative.le(t)).sum()),'denominator':int(eligible.sum()),'percent':100*relative[eligible].le(t).mean()} for t in [2,5,10]])
    csv(threshold,'figure4_nominal_difference_thresholds.csv')
    cov=[]
    independent=[(['homo_d','lumo_d','homo_a','lumo_a'],'Four energy levels'),(['d_a_ratio'],'D:A ratio'),(['solvent'],'Solvent'),(['device_structure'],'Device structure')]
    for fields,label in independent:
        available=b[fields].apply(base.present).all(axis=1)
        for key,m in masks.items():cov.append({'subset':key,'information':label,'fields':';'.join(fields),'records':int((m&available).sum()),'denominator':int(m.sum()),'percent':100*available[m].mean(),'used_as_selection_rule':False})
    cov=pd.DataFrame(cov);csv(cov,'figure4_independent_field_coverage.csv')
    rules['application_to_current_binary']={'input':'data/binary.csv','records':len(b),'purpose':'Recomputed analysis views; does not overwrite the current public data or claim row-level accuracy.',
        'nominal_pin_mw_cm2':100,'nominal_not_observed':'100 is a common normalization for this analysis, not a recovered per-record irradiance.',
        'ff_rescaled_records':int((complete&b.ff_num.le(1.5)).sum()),'complete_quartets':int(complete.sum()),'eligible_nominal_comparisons':int(eligible.sum()),
        'before_duplicate_rule':int(candidates.sum()),'removed_by_duplicate_rule':int((candidates&~strict).sum()),
        'duplicate_key_encoding':'Trimmed public metric strings and names, canonical fallback to reported names; numeric legacy IDs first in ascending order, followed by new nonnumeric IDs in lexical order. Keys do not include preparation conditions.',
        'missing_token_convention':sorted(base.EMPTY),'selection_assessment':'Compliance with required fields or the <=2% threshold is by construction, not independent validation.'}
    (SOURCE/'figure4_subset_rules.json').write_text(json.dumps(rules,ensure_ascii=False,indent=2),encoding='utf-8')
    return inputs,summary,threshold,cov

def load():
 b=base.read(DATA/'binary.csv');s=base.read(DATA/'binary_supplement.csv')
 b=b.merge(s.drop(columns='doi'),on='id',validate='one_to_one')
 b['doi_norm']=b['doi'].map(base.doi)
 meta=base.read(Path(__file__).parent/'assets/figure2_article_metadata.csv')
 b=b.merge(meta[['doi_norm','publish_year']],on='doi_norm',how='left',validate='many_to_one')
 b['publish_year']=pd.to_numeric(b.publish_year,errors='coerce')
 for f in ['pce','voc','jsc','ff']:b[f+'_num']=base.number(b[f])
 aliases={'PM6':['PM6','PBDB-TF','PBDB-T-2F'],'PTB7-Th':['PTB7-Th','PCE10','PCE-10','PBDTTT-EFT'],'PC71BM':['PC71BM','PC70BM','[70]PCBM'],'PC61BM':['PC61BM','PC60BM','[60]PCBM'],'IT-4F':['IT-4F','ITIC-4F'],'Y6':['Y6','BTP-4F','BTP-4F-8'],'Y6-BO':['Y6-BO','BTP-4F-12']}
 norm=lambda s:re.sub(r'\s+','',s).upper().replace('−','-').replace('–','-')
 mapping={norm(v):k for k,values in aliases.items() for v in values}
 alias_rows=[]
 for role in ['donor','acceptor']:
  b[role+'_label']=b[role].map(lambda s:mapping.get(norm(s),s.strip()))
  for r in b[['id','doi',role,role+'_label']].itertuples(index=False,name=None):
   if r[2].strip() and r[2].strip()!=r[3]:alias_rows.append({'id':r[0],'doi':r[1],'role':role,'reported_label':r[2],'plot_group':r[3]})
 csv(pd.DataFrame(alias_rows),'figure2_alias_groups.csv')
 return b

def fig2(b):
 dated=b.dropna(subset=['publish_year']).copy();dated.publish_year=dated.publish_year.astype(int)
 y=dated.groupby('publish_year').agg(records=('id','size'),source_dois=('doi_norm','nunique')).reset_index()
 g=b.groupby('doi_norm').size();bins=[]
 for label,lo,hi in [('1',1,1),('2',2,2),('3–5',3,5),('6–10',6,10),('11–20',11,20),('>20',21,np.inf)]:bins.append((label,int(g.between(lo,hi).sum())))
 csv(y,'figure2_yearly.csv')
 csv(g.rename('records').reset_index(),'figure2_records_per_doi.csv')
 csv(pd.DataFrame(bins,columns=['records_per_doi_bin','source_dois']),'figure2_record_bins.csv')
 material_rankings=[]
 fig,axes=plt.subplots(2,2,figsize=(7.09,5.6));fig.subplots_adjust(left=.15,right=.94,top=.92,bottom=.13,wspace=.65,hspace=.65)
 ax=axes[0,0];ax.bar(y.publish_year,y.records,color=C['blue']);ax.set_ylabel('Observation records');ax.set_xlabel('Publication year');ax.grid(axis='y')
 twin=ax.twinx();twin.plot(y.publish_year,y.source_dois,'o-',color=C['teal'],ms=3);twin.set_ylabel('Source DOIs',color=C['teal'],fontsize=8);twin.tick_params(axis='y',labelsize=7,colors=C['teal']);twin.yaxis.set_major_locator(MaxNLocator(4))
 ax.set_xticks([2011,2016,2021,2026]);ax.set_xticklabels(['2011','2016','2021','2026*']);style.panel(ax,'a','Literature coverage')
 ax=axes[0,1];ax.bar([v[0] for v in bins],[v[1] for v in bins],color=C['blue']);ax.set_ylim(0,3900);ax.set_ylabel('Source DOIs');ax.set_xlabel('Records per source DOI');ax.grid(axis='y')
 for i,(_,n) in enumerate(bins):ax.text(i,n+60,f'{n:,}',ha='center',fontsize=7.5)
 style.panel(ax,'b','Multiple observations per source');ax.text(.97,.96,f'{100*g.ge(2).mean():.1f}% have ≥2 records',transform=ax.transAxes,ha='right',va='top',fontsize=7.5)
 for ax,role,letter in zip(axes[1],['donor','acceptor'],'cd'):
  q=b[base.present(b[role+'_label'])].groupby(role+'_label').doi_norm.nunique().reset_index(name='source_dois').rename(columns={role+'_label':'label'})
  q=q[~q.label.map(lambda x:style.short_label(x) or x.strip().upper()=='PCBM')]
  material_rankings.append(q.assign(role=role))
  ranking(ax,q,'label',letter,'Frequent '+role+' labels',role=='acceptor')
  ax.tick_params(axis='y',labelsize=8)
 csv(pd.concat(material_rankings,ignore_index=True),'figure2_material_labels.csv')
 style.note(fig,'*2026 is partial. 72 records from 24 DOIs lack a year. Material rankings count source DOIs and omit ambiguous short labels.',x=.06,y=.025,width=108)
 style.save(fig,NAMES[1],'Figure 2. Literature and material coverage. a, Observation-record and source-DOI counts by publication year; 2026 is partial and 72 records from 24 DOIs lack a reliable year. b, Records per DOI; 91.3% of the 7,360 sources contribute at least two observations. c,d, Frequent donor and acceptor materials ranked by source-DOI count after harmonising recognised aliases and taking the DOI union for each material; ambiguous short labels and unspecified PCBM are omitted. Teal acceptor bars identify fullerene-name matches. Fullerene frequencies reflect historical studies and control devices; recent high-efficiency advances are led by non-fullerene acceptors.')

def fig3(b):
 perf=b[['pce_num','voc_num','jsc_num','ff_num']].notna().all(axis=1)
 present=lambda fs:b[fs].apply(base.present).all(axis=1)
 # Energy pairs require BOTH fields on the same observation; do not average coverage.
 fields=[('pce','PCE','Performance'),('voc',r'$V_{OC}$','Performance'),('jsc',r'$J_{SC}$','Performance'),('ff','FF','Performance'),
         ('donor','Donor name','Materials'),('acceptor','Acceptor name','Materials'),
         ('donor_smiles','Donor SMILES','Materials'),('acceptor_smiles','Acceptor SMILES','Materials'),
         ('homo_d+lumo_d','Donor HOMO/LUMO','Materials'),('homo_a+lumo_a','Acceptor HOMO/LUMO','Materials')]
 core_rows=[]
 for field,label,group in fields:
  mask=b[field+'_num'].notna() if group=='Performance' else present(field.split('+'))
  core_rows.append({'field':field,'label':label,'group':group,'records':int(mask.sum()),'denominator':len(b),'percent':100*mask.mean(),'coverage_rule':'finite numeric' if group=='Performance' else 'both nonempty' if '+' in field else 'nonempty'})
 core=pd.DataFrame(core_rows)
 core.to_csv(SOURCE/'figure3_core_field_coverage.csv',index=False,encoding='utf-8-sig',float_format='%.8g')
 device_fields=[('device_structure','Device structure','Device'),('device_type','Device type','Device'),('etl','ETL','Device'),('htl','HTL','Device'),('d_a_ratio','D:A ratio','Preparation'),('solvent','Solvent','Preparation')]
 device=pd.DataFrame([{'field':field,'label':label,'group':group,'records':int(present([field]).sum()),'denominator':len(b),'percent':100*present([field]).mean(),'coverage_rule':'nonempty'} for field,label,group in device_fields])
 device.to_csv(SOURCE/'figure3_device_preparation_coverage.csv',index=False,encoding='utf-8-sig',float_format='%.8g')
 names=present(['donor','acceptor']);sm=present(['donor_smiles','acceptor_smiles']);energy=present(['homo_d','lumo_d','homo_a','lumo_a'])
 specs=[('Four performance values',perf),('+ donor and acceptor names',perf&names),('+ paired SMILES',perf&sm),('+ donor and acceptor HOMO/LUMO',perf&energy)]
 rows=[{'combination':label.replace('\n',' '),'records':int(mask.sum()),'denominator':len(b),'percent':100*mask.mean()} for label,mask in specs]
 pd.DataFrame(rows).to_csv(SOURCE/'figure3_joint_availability.csv',index=False,encoding='utf-8-sig')
 # Selected preview 03: three field panels above four parallel joint-information plots.
 # The explanatory definitions belong in the caption, with no footer inside the figure.
 fig=plt.figure(figsize=(7.09,4.3))
 def axis(rect,label):
  ax=fig.add_axes(rect);ax.set_xlim(0,100)
  ax.set_xticks([0,25,50,75,100]);ax.set_xlabel(label,fontsize=7.2,labelpad=6)
  ax.spines['left'].set_visible(False);ax.spines['bottom'].set_color('#667681')
  ax.tick_params(axis='y',length=0,pad=6);ax.tick_params(axis='x',labelsize=7.1,pad=4)
  ax.grid(axis='x',color='#dce4e8',alpha=.85,lw=.5)
  return ax
 panels=[axis([x,.475,width,.385],'Records with information (%)')
         for x,width in [(.062,.205),(.445,.180),(.794,.178)]]
 def draw(ax,labels,values,color):
  yy=np.linspace(0,5,len(labels));ax.set_ylim(5.65,-.65)
  ax.barh(yy,values,height=.52,color=color,edgecolor='none',zorder=2)
  ax.set_yticks(yy);ax.set_yticklabels(labels,fontsize=7.5,linespacing=1.15)
  for y,value in zip(yy,values):
   ax.text(value+3,y,f'{value:.1f}',ha='left',va='center',fontsize=7.4,clip_on=False)
 performance=core[core.group.eq('Performance')];materials=core[core.group.eq('Materials')]
 draw(panels[0],performance.label,performance.percent,C['blue'])
 material_labels=[s.replace(' HOMO/LUMO','\nHOMO/LUMO') for s in materials.label]
 draw(panels[1],material_labels,materials.percent,C['teal'])
 draw(panels[2],device.label,device.percent,'#b0a99e')
 joint_labels=['Four performance\nvalues','+ donor/acceptor\nnames','+ paired SMILES','+ donor/acceptor\nHOMO/LUMO']
 for j,(label,row) in enumerate(zip(joint_labels,rows)):
  x=.026+.253*j;width=.198
  ax=axis([x,.15,width,.072],'Records with joint information (%)')
  ax.set_ylim(-.5,.5);ax.set_yticks([])
  ax.barh([0],[row['percent']],height=.40,color=C['blue'] if j==0 else C['teal'],edgecolor='none',zorder=2)
  ax.text(row['percent']+3,0,f"{row['percent']:.1f}",ha='left',va='center',fontsize=7.4,clip_on=False)
  fig.text(x+width/2,.253,label,ha='center',va='center',fontsize=7.6,linespacing=1.15)
  fig.text(x+width/2,.021,f"n = {row['records']:,}",ha='center',va='bottom',fontsize=7.5)
 for x,letter,title in zip([.022,.331,.687],'abc',['Performance','Materials','Device and preparation']):
  fig.text(x,.91,f'{letter}   {title}',fontsize=8.9,fontweight='bold')
 fig.add_artist(Line2D([.016,.984],[.355,.355],transform=fig.transFigure,color='#b9c6cd',lw=.6))
 fig.text(.022,.306,'d   Joint availability',fontsize=8.9,fontweight='bold')
 fig.text(.022,.974,'32,778 binary observation records',fontsize=8.7)
 style.save(fig,NAMES[2],'Figure 3. Availability of performance, material and device information. a, Core photovoltaic performance. b, Reported donor and acceptor names, mapped SMILES and material-reference energy values; donor and acceptor HOMO/LUMO each require both energy values for the same record. c, Device structure, device type, electron-transport layer (ETL), hole-transport layer (HTL), blend ratio and solvent. Device structure includes reported layer stacks or architecture descriptions. d, Four finite performance values alone, and combined separately with both reported material names, paired SMILES or all four donor/acceptor HOMO/LUMO values; n denotes the record count for each combination. Each combination is evaluated separately against the performance quartet. Every percentage uses all 32,778 binary records.')

def fig4(b):
 inputs,summary,threshold,cov=subset_analysis(b)
 fig,axes=plt.subplots(1,2,figsize=(7.09,3.5));fig.subplots_adjust(left=.10,right=.97,top=.77,bottom=.24,wspace=.65)
 ax=axes[0];im=ax.hexbin(inputs.pce_num,inputs.pce_calculated_nominal_percent,gridsize=40,mincnt=1,cmap=DENSITY_CMAP,norm=LogNorm(),linewidths=0)
 vmax=max(inputs.pce_num.max(),inputs.pce_calculated_nominal_percent.max())*1.05;ax.plot([0,vmax],[0,vmax],'--',lw=.7,color=C['ink']);ax.set_xlim(0,vmax);ax.set_ylim(0,vmax);ax.set_aspect('equal');ax.set_xlabel('Reported PCE (%)');ax.set_ylabel('Calculated PCE (%)');vector_colorbar(fig,im,ax,'Records / hexagon');style.panel(ax,'a','Full binary collection',pad=19);ax.text(0,1.03,f'n = {len(inputs):,} complete quartets',transform=ax.transAxes,fontsize=7.5)
 ax=axes[1];xmax=max(inputs.relative_difference_percent.max()*1.15,12)
 csv(pd.DataFrame({'reported_pce_center':im.get_offsets()[:,0],'calculated_pce_center':im.get_offsets()[:,1],'records':im.get_array()}),'figure4_nominal_parity_hexbin.csv')
 cumulative=[]
 for key,label,color,ls in [('full','Full binary',C['gray'],'-'),('strict_performance','Strict performance',C['blue'],'-'),('strict_molecular','Strict molecular',C['teal'],'--')]:
  q=inputs.loc[inputs[key]].sort_values(['relative_difference_percent','id'],kind='stable');xx=q.relative_difference_percent.to_numpy();yy=100*np.arange(1,len(q)+1)/len(q)
  cumulative.append(q[['id','relative_difference_percent']].assign(subset=key,denominator=len(q),cumulative_percent=yy))
  ax.step(np.r_[0,xx,xmax],np.r_[0,yy,100],where='post',color=color,ls=ls,lw=1.4,label=f'{label} (n = {len(q):,})')
 csv(pd.concat(cumulative,ignore_index=True),'figure4_nominal_cumulative_difference.csv')
 ax.set_xscale('symlog',linthresh=.1,linscale=.5);ax.set_xlim(0,xmax);ax.set_ylim(0,102);ticks=[x for x in [0,.1,1,10,100,1000] if x<=xmax];ax.set_xticks(ticks);ax.set_xticklabels([f'{x:g}' for x in ticks]);ax.set_yticks([0,25,50,75,100]);ax.grid();ax.set_xlabel('Absolute relative difference (%)');ax.set_ylabel('Cumulative records (%)')
 for threshold_value in [2,5,10]:ax.axvline(threshold_value,ls=':',lw=.6,color=C['gray'])
 ax.legend(fontsize=6.8,frameon=True,facecolor='white',edgecolor='none',loc='lower right');style.panel(ax,'b','Full data and selected subsets',pad=19)
 style.note(fig,'Nominal irradiance: 100 mW cm$^{-2}$. Strict subsets apply the 2% threshold; their curves describe selection, not independent validation.',x=.055,y=.025,width=112)
 style.save(fig,NAMES[3],'Figure 4. Photovoltaic consistency and subset selection. a, Reported versus calculated PCE for 32,254 complete positive-PCE binary observations, using nominal irradiance of 100 mW cm⁻²; hexagon colour shows logarithmic record counts and the dashed line indicates equality. b, Cumulative absolute relative differences for the full evaluable collection, strict performance subset and strict molecular subset, each normalised by its own size. The axis is linear below 0.1% and logarithmic above. Dotted vertical lines mark 2%, 5% and 10%. The strict subsets incorporate a 2% relative-difference threshold, and their curves show the resulting distributions.')

def render_workflow():
    from svglib.svglib import svg2rlg, register_font as register_svg_font
    from reportlab.graphics import renderPDF
    from reportlab.pdfbase import pdfmetrics
    from matplotlib.font_manager import FontProperties, findfont
    import fitz
    # Map the SVG's Arial family explicitly; Matplotlib supplies a portable
    # fallback when Arial is absent. No font files are redistributed.
    for name,weight,fontstyle in [('Arial','normal','normal'),('Arial-Bold','bold','normal'),('Arial-Italic','normal','italic')]:
        font=findfont(FontProperties(family='Arial',weight=weight,style=fontstyle))
        register_svg_font('Arial',font,weight=weight,style=fontstyle,rlgFontName=name)
    regular=findfont(FontProperties(family='Arial'))
    register_svg_font('Arial',regular,weight='400',rlgFontName='Arial')
    register_svg_font('sans-serif',regular,rlgFontName='Arial')
    pdfmetrics.registerFontFamily('Arial',normal='Arial',bold='Arial-Bold',italic='Arial-Italic',boldItalic='Arial-Bold')
    source=Path(__file__).parent/'assets/fig1_workflow_and_scope.svg'
    target=FIG/'fig1_workflow_and_scope.svg'
    shutil.copy2(source,target)
    drawing=svg2rlg(str(target))
    renderPDF.drawToFile(drawing,str(target.with_suffix('.pdf')))
    with fitz.open(target.with_suffix('.pdf')) as doc:
        doc[0].get_pixmap(dpi=300).save(str(target.with_suffix('.png')))

def validate():
    inventory=[]
    for table,expected_rows,expected_cols,expected_dois in [('binary',32778,37,7360),('binary_supplement',32778,47,7360),('ternary',2588,68,461),('indoor',1678,57,166),('materials_reference',4988,7,None),('binary_strict_performance',27863,37,6873),('binary_strict_molecular',19299,37,5193)]:
        df=base.read(DATA/(table+'.csv'))
        assert df.shape==(expected_rows,expected_cols),(table,df.shape)
        if expected_dois is not None:
            assert base.present(df.doi).all()
            assert df.doi.map(base.doi).nunique()==expected_dois
        inventory.append({'table':table,'rows':len(df),'columns':len(df.columns),'source_dois':expected_dois})
    b=base.read(DATA/'binary.csv');s=base.read(DATA/'binary_supplement.csv')
    assert b.id.is_unique and s.id.is_unique
    assert b.set_index('id').doi.sort_index().equals(s.set_index('id').doi.sort_index())
    header=json.loads((Path(__file__).parent/'binary_original_columns.json').read_text())
    assert b.columns.tolist()==header,'The original 37-column schema must be preserved'
    csv(pd.DataFrame(inventory),'dataset_inventory.csv')

def supplement(data):
    t=data['ternary'];i=data['indoor']
    st=[(('component_1_name','component_2_name','component_3_name'),'Three names','text'),(('component_1_role','component_2_role','component_3_role'),'Three roles','text'),(('ratio_raw',),'Composition ratio','text'),(('pce_percent',),'PCE','number'),(('voc_v','jsc_ma_cm2','ff_percent','pce_percent'),'Four metrics','number')]
    si=[(('light_source_type',),'Light source','text'),(('illuminance_lux',),'Illuminance (lux)','number'),(('irradiance_uw_cm2',),'Irradiance','number'),(('pce_percent',),'PCE','number'),(('voc_v',),r'$V_{OC}$','number'),(('jsc_ua_cm2',),r'$J_{SC}$','number'),(('ff_percent',),'FF','number'),(('voc_v','jsc_ua_cm2','ff_percent','pce_percent'),'Four metrics','number')]
    ct=style.coverage(t,st,'ternary');ci=style.coverage(i,si,'indoor');csv(pd.concat([ct,ci]),'figureS1_collection_field_coverage.csv')
    fig,axes=plt.subplots(1,2,figsize=(style.WIDTH,3.3));fig.subplots_adjust(left=.18,right=.97,top=.78,bottom=.23,wspace=.83)
    for a,cov,color,letter,title in zip(axes,[ct,ci],[C['teal'],C['blue']],['a','b'],['Ternary collection\n(n = 2,588)','Indoor collection\n(n = 1,678)']):
        style.coverage_bars(a,cov,color);a.tick_params(labelsize=7);a.set_xlabel('Records with information (%)',fontsize=7.4);style.panel(a,letter,title,pad=12)
    style.save(fig,NAMES[4],'Supplementary Figure S1. Information available in the ternary and indoor collections. a, Three component names, three roles, reported composition ratio, PCE and complete numerical performance quartets among 2,588 ternary observations. b, Light-source identity, illuminance, irradiance and performance information among 1,678 indoor observations. Bars show percentages; adjacent numbers give record counts. Text fields use the documented missing-token convention; performance and illumination values require finite numbers. Illuminance and irradiance are reported separately. Both collections have partial literature coverage.')

def main():
    global RELEASE,DATA,FIG,SOURCE
    ap=argparse.ArgumentParser(description='Reproduce OPV-DB v2.0.0 statistics and main figures offline.')
    ap.add_argument('--release',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();RELEASE=args.release.resolve();DATA=RELEASE/'data'
    output=args.out.resolve()
    if output==RELEASE or output in [DATA,RELEASE/'code',RELEASE/'source_data']:
        ap.error('Choose a separate output directory.')
    FIG=output/'figures';SOURCE=output/'source_data'
    FIG.mkdir(parents=True,exist_ok=True);SOURCE.mkdir(parents=True,exist_ok=True)
    style.FIG=FIG;style.C=C;style.setup();validate();b=load()
    render_workflow();fig2(b);fig3(b);fig4(b)
    supplement({name:base.read(DATA/(name+".csv")) for name in ["ternary","indoor"]})
    joint=base.read(SOURCE/'figure3_joint_availability.csv');assert joint.records.astype(int).tolist()==[32254,31677,22287,20973]
    subsets=base.read(SOURCE/'figure4_subset_summary.csv');assert subsets.records.astype(int).tolist()==[32778,27863,19299]
    membership=base.read(SOURCE/'figure4_current_subset_membership.csv')
    raw=base.read(DATA/'binary.csv')
    for key in ['strict_performance','strict_molecular']:
        expected=raw[membership[key].str.lower().eq('true').to_numpy()].reset_index(drop=True)
        supplied=base.read(DATA/('binary_'+key+'.csv')).reset_index(drop=True)
        assert supplied.equals(expected), 'Published subset differs from recomputed membership/cells: '+key
    (output/'figure_captions.json').write_text(json.dumps(style.CAPTIONS,ensure_ascii=False,indent=2),encoding='utf-8')
    result={'status':'passed','binary_rows':len(b),'binary_columns':37,'binary_source_dois':b.doi_norm.nunique(),'joint_availability_records':joint.records.astype(int).tolist(),'subset_records':subsets.records.astype(int).tolist(),'network_required':False}
    (output/'verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result))
if __name__=='__main__':main()
