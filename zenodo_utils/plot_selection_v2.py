#!/usr/bin/env python3
"""Plot v2 converter sidecars; no ROOT or shower-tensor rescan needed.

CLD compares the v1 part1-4 RULE (inclusive), v1 part5-10 RULE (numMC==1),
and v2 on exactly the SAME input ROOT pool. These are policy comparisons,
not claims that all historical partitions were checked or regenerated.
For ALLEGRO labels are inclusive reference / v1 numMC==1 / v2.
Each chart is generated separately; overview PNGs are contact sheets of these
charts. Endpoint plots show primary ENDPOINTS, never secondary vertices.
"""
from __future__ import annotations
import argparse
import csv
import html
import json
import math
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from PIL import Image
# Use the converter from this checkout, independently of the working directory.
# Importing it does not import ROOT or start conversion.
_CONVERTER_DIR = Path(__file__).resolve().parents[1] / 'scripts'
if not (_CONVERTER_DIR / 'edm4hep2h5_v2.py').is_file():
    raise SystemExit('Missing scripts/edm4hep2h5_v2.py in this ddfastsim checkout')
sys.path.insert(0, str(_CONVERTER_DIR))
from edm4hep2h5_v2 import GEOMETRIES, support


KEYS = ('incident_energy_MeV','incident_theta','incident_phi','etot_MeV',
        'endpoint_x_mm','endpoint_y_mm','endpoint_z_mm','nmc','reason')
CHARTS = ('01_energy_spectrum','02_costheta_spectrum','03_phi_spectrum',
          '04_total_energy_spectrum','05_response','06_response_vs_energy',
          '07_primary_endpoints','08_efficiency_vs_energy')


def read_inputs(paths):
    blocks={k:[] for k in KEYS}; sources=[]; detector=None; seen=set(); settings=None
    for path in paths:
        with np.load(path,allow_pickle=False) as f:
            meta=json.loads(str(f['metadata_json'].item()))
            if meta.get('status')!='complete': raise ValueError('Incomplete conversion: '+str(path))
            conf=meta['settings']
            if detector is None:
                detector=conf['detector']; settings=conf
            if conf != settings:
                raise ValueError('Do not combine different detector/settings/converter versions')
            identity=meta['input']['path']
            if identity in seen: raise ValueError('Duplicate ROOT input: '+identity)
            seen.add(identity); n=meta['total_events']
            if any(f[k].shape != (n,) for k in KEYS): raise ValueError('Malformed cache '+str(path))
            for k in KEYS: blocks[k].append(f[k].copy())
            sources.append(meta)
    if not sources: raise ValueError('No selection sidecars found')
    return {k:np.concatenate(v) for k,v in blocks.items()},detector,settings,sources


def make_bins(values,n=60,lower=None,upper=None):
    finite=values[np.isfinite(values)]
    if finite.size==0: raise ValueError('Cannot bin empty/nonfinite input')
    lo=float(finite.min()) if lower is None else lower
    hi=float(finite.max()) if upper is None else upper
    if hi<=lo: lo-=max(abs(lo)*.01,.001); hi+=max(abs(hi)*.01,.001)
    return np.linspace(lo,np.nextafter(hi,np.inf),n+1)


def hist_density(v,bins):
    finite=np.isfinite(v); count,_=np.histogram(v[finite],bins)
    # Divisor includes all finite records, not just the ones in a displayed range.
    return count/(max(1,int(finite.sum()))*np.diff(bins))


def binned_mean(x,y,mask,bins):
    n=np.zeros(len(bins)-1,int); mean=np.full(len(n),np.nan); sem=mean.copy()
    b=np.searchsorted(bins,x,side='right')-1
    for i in range(len(n)):
        values=y[mask & (b==i) & np.isfinite(y)]
        n[i]=len(values)
        if len(values): mean[i]=values.mean()
        if len(values)>1: sem[i]=values.std(ddof=1)/math.sqrt(len(values))
    return n,mean,sem


def lines(ax,g):
    # Nominal REGION envelopes, not just active absorber/sensitive material.
    for x1,x2,z1,z2,label in (
        (g['barrel_r'],g['barrel_outer_r'],0,g['barrel_half_z'],'Barrel region'),
        (g['endcap_inner_r'],g['endcap_outer_r'],g['endcap_z'],g['endcap_back_z'],'Endcap region')):
        ax.plot([x1,x2,x2,x1,x1],[z1,z1,z2,z2,z1],linestyle='--',linewidth=1.2,label=label)
    ax.axhline(2000,linestyle=':',linewidth=1.,label='|z|=2000 mm analysis divider')
    ax.set_xlabel('Primary endpoint polygon support (mm)' if g['sides'] else 'Primary endpoint radius (mm)')
    ax.set_ylabel(r'Primary endpoint $|z|$ (mm)')


def save(fig,path):
    fig.tight_layout()
    fig.savefig(path,dpi=160)
    plt.close(fig)


def contact_sheet(paths,out):
    """Assemble already-generated separate charts; never resample event populations."""
    images=[Image.open(p).convert('RGB') for p in paths]
    w,h=images[0].size
    # Use the source figure background rather than specifying a chart colour.
    sheet=Image.new('RGB',(2*w,4*h),images[0].getpixel((0,0)))
    for i,img in enumerate(images):
        sheet.paste(img,(i%2*w,i//2*h))
        img.close()
    sheet.save(out)


def run(args):
    paths=[]
    for p in args.inputs:
        p=Path(p)
        paths.extend(sorted(p.glob('*.selection.npz')) if p.is_dir() else [p])
    d,detector,settings,sources=read_inputs(paths)
    g=settings['geometry']
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists() and not args.replace:
        raise FileExistsError('Plot output exists; use --replace to regenerate plots only')
    E=d['incident_energy_MeV'].astype(float)/1000.
    T=d['etot_MeV']/1000.
    R=np.divide(T,E,out=np.full(len(E),np.nan),where=E>0)
    valid=np.isfinite(E)&(E>0)&np.isfinite(d['incident_theta'])&np.isfinite(d['incident_phi'])
    if not valid.all():
        raise ValueError('Invalid conditioning found; inspect conversion reports before comparing spectra')
    c=np.cos(d['incident_theta'].astype(float)); phi=d['incident_phi'].astype(float)
    initial=np.ones(len(E),bool)
    labels=(['v1 part1-4 rule: no cut','v1 part5-10 rule: numMC==1','v2: geometry + early-time veto']
            if detector=='CLD' else ['Unfiltered reference','v1 rule: numMC==1','v2: geometry + early-time veto'])
    names=('inclusive','numMC1','v2')
    masks=(initial,d['nmc']==1,d['reason']==0)
    Eb=make_bins(E,20); Ehist=make_bins(E)
    cb=make_bins(c); pb=np.linspace(-np.pi,np.nextafter(np.pi,np.inf),61)
    Tb=make_bins(T,lower=0.); Rb=make_bins(R,80,lower=0.)
    face=np.array([support(x,y,g) for x,y in zip(d['endpoint_x_mm'],d['endpoint_y_mm'])])
    absz=np.abs(d['endpoint_z_mm'])
    finite_end=np.isfinite(face)&np.isfinite(absz)
    xmax=max(g['barrel_outer_r'],g['endcap_outer_r'],float(np.nanmax(face)))*1.025
    zmax=max(g['endcap_back_z'],float(np.nanmax(absz)))*1.025
    endpoint_bins=(np.linspace(0,xmax,161),np.linspace(0,zmax,161))
    allcounts,_,_=np.histogram2d(face[finite_end],absz[finite_end],bins=endpoint_bins)
    vmax=max(2.,float(allcounts.max()))
    metrics=[(CHARTS[0],E,Ehist,r'$E_{\rm gen}$ (GeV)'),
             (CHARTS[1],c,cb,r'$\cos\theta_{\rm gen}$'),
             (CHARTS[2],phi,pb,r'$\phi_{\rm gen}$ (rad)'),
             (CHARTS[3],T,Tb,r'$E_{\rm tot}$ in stored voxels (GeV)'),
             (CHARTS[4],R,Rb,r'$E_{\rm tot}/E_{\rm gen}$')]
    counts0,_=np.histogram(E,Eb)
    overview_paths=[]; summary={}; csv_rows=[]
    policy_info=list(zip(names,labels,masks))
    # Separate chart files plus one overview contact sheet per rule and combined.
    for panel in names+('comparison',):
        pd=out/panel; pd.mkdir(exist_ok=True)
        active=policy_info if panel=='comparison' else [p for p in policy_info if p[0]==panel]
        title=detector+' | '+('Selection-rule comparison' if panel=='comparison' else active[0][1])
        for filename,values,bins,xlabel in metrics:
            fig,ax=plt.subplots(figsize=(7.4,5.2))
            for name,label,mask in active:
                ax.stairs(hist_density(values[mask],bins),bins,label=f'{label} (N={mask.sum():,})')
            ax.set(xlabel=xlabel,ylabel='Probability density',title=title)
            ax.legend(fontsize=8)
            save(fig,pd/(filename+'.png'))
        fig,ax=plt.subplots(figsize=(7.4,5.2))
        for name,label,mask in active:
            nn,mm,ss=binned_mean(E,R,mask,Eb)
            ax.errorbar((Eb[:-1]+Eb[1:])/2,mm,yerr=ss,marker='.',linestyle='-',label=label)
        ax.set(xlabel=r'$E_{\rm gen}$ (GeV)',ylabel=r'Mean $E_{\rm tot}/E_{\rm gen}$',title=title)
        ax.legend(fontsize=8); save(fig,pd/(CHARTS[5]+'.png'))
        fig,ax=plt.subplots(figsize=(7.4,5.2))
        if panel!='comparison':
            m=active[0][2]&finite_end
            counts,xx,zz=np.histogram2d(face[m],absz[m],bins=endpoint_bins)
            im=ax.pcolormesh(xx,zz,np.ma.masked_less(counts.T,1),norm=LogNorm(vmin=1,vmax=vmax),rasterized=True)
            fig.colorbar(im,ax=ax,label='Events per bin (common colour scale)')
        else:
            # All events, no unequal downsampling. Nested populations overlap.
            for (name,label,mask),marker,sz in zip(active,('x','+','.'),(8,5,1.5)):
                m=mask&finite_end
                ax.scatter(face[m],absz[m],s=sz,marker=marker,alpha=.35,label=label,rasterized=True)
        lines(ax,g); ax.set(xlim=(0,xmax),ylim=(0,zmax),title=title)
        ax.legend(fontsize=7,loc='upper right')
        save(fig,pd/(CHARTS[6]+'.png'))
        fig,ax=plt.subplots(figsize=(7.4,5.2))
        for name,label,mask in active:
            nn,_=np.histogram(E[mask],Eb)
            eff=np.divide(nn,counts0,out=np.full(len(nn),np.nan),where=counts0>0)
            ax.stairs(eff,Eb,label=label)
        ax.set(xlabel=r'$E_{\rm gen}$ (GeV)',ylabel='Selected / original events in bin',ylim=(0,1.025),title=title)
        ax.legend(fontsize=8);save(fig,pd/(CHARTS[7]+'.png'))
        overview=out/(panel+'_overview.png')
        contact_sheet([pd/(key+'.png') for key in CHARTS],overview)
        overview_paths.append((title,overview.name,panel))
    for name,label,mask in policy_info:
        vals=R[mask]; finite=vals[np.isfinite(vals)]
        summary[name]=dict(label=label,kept=int(mask.sum()),fraction=float(mask.mean()),
                           mean_Egen_GeV=float(E[mask].mean()) if mask.any() else None,
                           mean_Etot_GeV=float(T[mask].mean()) if mask.any() else None,
                           mean_response=float(finite.mean()) if len(finite) else None,
                           response_quantiles_05_50_95=np.quantile(finite,[.05,.5,.95]).tolist() if len(finite) else [],
                           zero_grid_events=int(np.sum(mask&(T==0))))
        for variable,values,bins in [('energy_GeV',E,Eb),('costheta',c,make_bins(c,20)),('phi_rad',phi,np.linspace(-np.pi,np.nextafter(np.pi,np.inf),21))]:
            base,_=np.histogram(values,bins); kept,_=np.histogram(values[mask],bins)
            for i in range(len(base)):
                csv_rows.append([name,variable,bins[i],bins[i+1],int(base[i]),int(kept[i]),float(kept[i]/base[i]) if base[i] else ''])
    output=dict(detector=detector,total_ROOT_events=len(E),files=[s['input']['path'] for s in sources],
                settings=settings,scenarios=summary,
                meaning='Same original ROOT population, three selection rules. CLD v1 labels identify historic RULES, not measurements of independent historical partition sets.',
                endpoints='Original photon endpoints, never secondary vertices; every event contributes. Dashed envelopes are nominal region boundaries.',
                caveats=['Endpoint predicate does not establish unchanged ECAL entry state.',
                         'Response means are at each selection\'s own input distribution; curves vs energy use finite bins.',
                         'Joint endpoint clouds are nested and overlap; individual density plots share bins and colour normalization.',
                         'No claims of duplicate-free or statistically independent events.',
                         'Uniformity plots do not correct for selection; cos(theta), not theta, is plotted.'])
    (out/'summary.json').write_text(json.dumps(output,indent=2)+'\n')
    with open(out/'selection_efficiencies.csv','w',newline='') as f:
        w=csv.writer(f);w.writerow(['scenario','variable','bin_low','bin_high','original','selected','efficiency']);w.writerows(csv_rows)
    chunks=['<!doctype html><meta charset="utf-8"><title>LEMURS v2 selection</title>',
            '<style>body{font:16px sans-serif;max-width:1500px;margin:2em auto}img{max-width:100%}section{margin:3em 0}</style>',
            f'<h1>{detector}: three selection rules on {len(E):,} original ROOT events</h1>',
            '<p>CLD v1 labels describe the historic selection rules, applied here to the same input events. '
            'ALLEGRO uses its own region dimensions. No events are regenerated. '
            'Dashed lines are nominal ECAL region envelopes; the dotted |z|=2000 mm line is only an analysis divider.</p>']
    for title,png,panel in overview_paths:
        chunks.append(f'<section><h2>{html.escape(title)}</h2><a href="{png}"><img src="{png}"></a></section>')
    (out/'index.html').write_text('\n'.join(chunks))
    print(json.dumps({'detector':detector,'total_ROOT_events':len(E),'scenarios':summary,'plots':str(out.resolve())},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs',nargs='+',help='*.selection.npz files or a detector output directory')
    p.add_argument('--outdir',required=True)
    p.add_argument('--replace',action='store_true',help='Replace previously generated plots, never data')
    try: run(p.parse_args())
    except (ValueError,OSError,KeyError) as e:
        print('ERROR: '+str(e),file=sys.stderr);sys.exit(1)
