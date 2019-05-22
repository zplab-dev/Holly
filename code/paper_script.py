from elegant import worm_data
from elegant import load_data
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from zplib.image import colorize
import numpy
import scipy.stats
import time
import statsmodels.api as sm
from pandas import DataFrame
from sklearn import linear_model
from glob import glob
import re
import pathlib
from zplib import datafile
import freeimage
from elegant import worm_spline
import os

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['figure.figsize'] = [15.0,9.0]
plt.rcParams['xtick.direction']='inout'
plt.rcParams['ytick.direction']='inout'
plt.rcParams['axes.spines.top']=False
plt.rcParams['axes.spines.right']=False
plt.rcParams['lines.dash_capstyle']='round'
plt.rcParams['lines.solid_capstyle']='round'
plt.rcParams['savefig.transparent']=True
plt.rcParams['legend.frameon']=False

plt.rcParams['axes.labelpad']=15.0
plt.rcParams['savefig.transparent']=True
plt.rcParams['legend.labelspacing']=.5
plt.rcParams['lines.linewidth']=3
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rcParams['xtick.major.pad'] = 9
plt.rcParams['ytick.major.pad'] = 9
plt.rcParams['lines.markersize']= 10

TIME_UNITS = dict(days=24, hours=1, minutes=1/60, seconds=1/3600)

label_fontdict={'size':20,'family':'arial'}
title_fontdict={'size':22, 'weight':'bold','family':'arial'}
text_fontdict={'size':20,'weight':'bold','family':'arial'}
gfp_measures={
'gfp_mean':'mean intensity',
'gfp_median':'median intensity',
'gfp_maximum':'maximum intensity',
'gfp_over_99_mean':'mean of pixels over 99th percentile intensity',
'gfp_percentile_95':'95th percentile intensity',
'gfp_percentile_99':'99th percentile intensity',
'gfp_sum':'summed intensity',
'green_yellow_excitation_autofluorescence_percentile_95': '95th percentile intensity of autofluorescence'}


#generates single trace of average expression for each measure in gfp_measures
def figure_1_population_trajectory(root_folder,overwrite=True,min_age=24,time_units='hours'):
	subdirectories=glob(root_folder+'/*')
	for subdirectory in subdirectories:

		exp_dirs=glob(subdirectory+'/*/')
		worms=process_worms(exp_dirs, prefixes=[dir+' ' for dir in exp_dirs])

		if pathlib.Path(subdirectory+'/'+'parameters.tsv').exists():
			header,parameters=datafile.read_delimited(subdirectory+'/'+'parameters.tsv')
			parameters=list(parameters)
			max_age=parameters[header.index('max_age')]
			max_age=((int)(max_age[0]))*24
		else:
			max_age=numpy.inf		

		for key,value in gfp_measures.items():
			miRNA=re.findall(r"\w\w\w-\d+",subdirectory)[0]
			save_name='/Volumes/9karray/Kinser_Holly/all_figures_figure1/'+miRNA+'_'+key+'.svg'
			if (pathlib.Path(save_name).exists()==True) and (overwrite==False):
				continue
			try:
				worms.get_time_range(key,24,numpy.inf)
			except:
				continue	
			
			if len(exp_dirs)>1:
				median_gfp={}
				groups=worms.group_by([w.name.split()[0] for w in worms])	
	
				for exp, wormies in groups.items():
						data = wormies.get_time_range(key,min_age,max_age)	

						median_gfp[exp] = numpy.median([d[:,1].mean() for d in data])
				for worm in worms:
						exp = worm.name.split()[0]
						worm.td.scaled_gfp = getattr(worm.td,key)/ median_gfp[exp] * median_gfp[exp_dirs[0]]
						

			else:
				for worm in worms:

						worm.td.scaled_gfp = getattr(worm.td,key)

			lifespans=worms.get_feature('lifespan')
			bins={'['+(str)(lifespans.min())+'-'+(str)(lifespans.max())+']':worms}



			averaged_worms = worm_data.meta_worms(bins, 'scaled_gfp', age_feature='age')
			figure=averaged_worms.plot_timecourse('scaled_gfp',age_feature='age',min_age=min_age,max_age=max_age,time_units=time_units)
			
			cohorts=[]

			for mew, mewmew in bins.items():
				cohorts.append(len(mewmew))	
			
			text='n = ' + (str)(cohorts)
			plt.figtext(.5, .15, text, fontdict=text_fontdict)
			plt.title('Average ' +'P'+miRNA+'::GFP expression profiles',fontdict=title_fontdict)
			plt.ylabel('P'+miRNA+'::GFP expression ('+value+')',fontdict=label_fontdict)
			print('Saving '+save_name)	
			plt.savefig(save_name)
			plt.gcf().clf
			plt.close()
#generates cohort traces of average expression for each measure in gfp_measures
def figure_2_population_trajectory_cohorts(root_folder,overwrite=True,min_age=24,time_units='hours'):
	subdirectories=glob(root_folder+'/*')
	for subdirectory in subdirectories:

		exp_dirs=glob(subdirectory+'/*')
		worms=process_worms(exp_dirs, prefixes=[dir+' ' for dir in exp_dirs])

		if pathlib.Path(subdirectory+'/'+'parameters.tsv').exists():
			header,parameters=datafile.read_delimited(subdirectory+'/'+'parameters.tsv')
			parameters=list(parameters)
			max_age=parameters[header.index('max_age')]
			max_age=((int)(max_age[0]))*24
		else:
			max_age=numpy.inf	
		
		for key,value in gfp_measures.items():
			miRNA=re.findall(r"\w\w\w-\d+",subdirectory)[0]
			save_name='/Volumes/9karray/Kinser_Holly/all_figures_figure2/'+miRNA+'_'+key+'.svg'
			if (pathlib.Path(save_name).exists()==True) and (overwrite==False):
				continue
			try:
				worms.get_time_range(key,min_age,numpy.inf)
			except:
				continue	
			
			if len(exp_dirs)>1:
				median_gfp={}
				groups=worms.group_by([w.name.split()[0] for w in worms])	
	
				for exp, wormies in groups.items():
						data = wormies.get_time_range(key,min_age,numpy.inf)	

						median_gfp[exp] = numpy.median([d[:,1].mean() for d in data])
				for worm in worms:
						exp = worm.name.split()[0]
						worm.td.scaled_gfp = getattr(worm.td,key)/ median_gfp[exp] * median_gfp[exp_dirs[0]]
						

			else:
				for worm in worms:

						worm.td.scaled_gfp = getattr(worm.td,key)

			lifespans=worms.get_feature('lifespan')
			bins=worms.bin('lifespan',nbins=4, equal_count=True)



			averaged_worms = worm_data.meta_worms(bins, 'scaled_gfp', age_feature='age')
			figure=averaged_worms.plot_timecourse('scaled_gfp',age_feature='age',min_age=min_age, max_age=max_age,time_units=time_units)
			
			cohorts=[]

			for mew, mewmew in bins.items():
				cohorts.append(len(mewmew))	
			
			text='n = ' + (str)(cohorts)

			plt.figtext(.5, .15, text, fontdict=text_fontdict)
			plt.title('Average ' +'P'+miRNA+'::GFP expression profiles',fontdict=title_fontdict)
			plt.ylabel('P'+miRNA+'::GFP expression ('+value+')',fontdict=label_fontdict)
			print('Saving '+save_name)	
			plt.savefig(save_name)
			plt.gcf().clf
			plt.close()
#regresses on slope and mean for each measure in gfp_measures, from 3 dph to timepoint at which 90% of population is alive
def figure_2_regression_on_slope_and_mean(root_folder,overwrite=True,min_age=3*24,target='lifespan'):
	subdirectories=glob(root_folder+'/*')
	for subdirectory in subdirectories:

		exp_dirs=glob(subdirectory+'/*')
		worms=process_worms(exp_dirs, prefixes=[dir+' ' for dir in exp_dirs])
		if pathlib.Path(subdirectory+'/'+'parameters.tsv').exists():
			header,parameters=datafile.read_delimited(subdirectory+'/'+'parameters.tsv')
			parameters=list(parameters)
			max_age=parameters[header.index('max_age')]
			max_age=((int)(max_age[0]))*24
		else:
			lifespans=worms.get_feature(target)
			max_age=round(numpy.percentile(lifespans,10))
		
		for key,value in gfp_measures.items():
			miRNA=re.findall(r"\w\w\w-\d+",subdirectory)[0]
			save_name='/Volumes/9karray/Kinser_Holly/all_figures_figure2b/'+miRNA+'_'+key+'.svg'
			
			if (pathlib.Path(save_name).exists()==True) and (overwrite==False):
				continue
		
			try:
				results=worms.regress(get_average_wrapper(min_age,max_age,key),get_slope_wrapper(min_age,max_age,key),target=target)
			except:	
				continue	
			
			color_vals = colorize.scale(lifespans, output_max=1)
			colors = colorize.color_map(color_vals, uint8=False)	
			plt.scatter(lifespans,results.y_est,c=colors)
			(m,b) = numpy.polyfit(lifespans, results.y_est, 1)
			yp=numpy.polyval([m,b], lifespans)
			mew, mewmew=zip(*sorted(zip(lifespans, yp)))
			mew, mewmew = (list(t) for t in zip(*sorted(zip(mew, mewmew))))
			plt.plot([mew[0],mew[-1]],[mewmew[0],mewmew[-1]],c='gray',alpha=.7)
			ftext="$r^{2}$ = "+(str)(round(results.R2,3))
			more_text="n= " + (str)(len(lifespans))
			plt.title('Regression on Slope and Mean ' +miRNA+'::GFP expression profiles '+(str)(min_age) + ' to '+(str)(max_age)+ ' '+key,fontdict=title_fontdict)
			plt.xlabel('Actual lifespan (days)',fontdict=label_fontdict)
			plt.ylabel("Predicted lifespan (days)",fontdict=label_fontdict)
			plt.figtext(.15,.8,ftext,fontsize=20,ha='left')
			plt.figtext(.8, .2, more_text, fontsize=20, ha='left')
			print('Saving '+save_name)	
			plt.savefig(save_name)
			plt.gcf().clf
			plt.close()
#regresses on slope and mean for each measure in gfp_measures, from 3 dph to timepoint at which 90% of population is alive. 
#Controls for autofluorescence in regression which eliminates spurious correlation with dim reporters			
def figure_2_regression_on_slope_and_mean_controlled_for_autofluorescence(root_folder,overwrite=True,min_age=3*24,target='lifespan'):
	subdirectories=glob(root_folder+'/*')
	for subdirectory in subdirectories:

		exp_dirs=glob(subdirectory+'/*')
		worms=process_worms(exp_dirs, prefixes=[dir+' ' for dir in exp_dirs])
		for key,value in gfp_measures.items():
			miRNA=re.findall(r"\w\w\w-\d+",subdirectory)[0]
			save_name='/Volumes/9karray/Kinser_Holly/all_figures_figure2b/'+miRNA+'_'+key+'.svg'
			if (pathlib.Path(save_name).exists()==True) and (overwrite==False):
				continue
		
			lifespans=worms.get_feature(target)
			max_age=round(numpy.percentile(lifespans,10))
			try:
				results=worms.regress(get_average_wrapper(min_age,max_age,key),get_slope_wrapper(min_age,max_age,key),target=target,control_features=[get_average_wrapper(min_age,max_age,'green_yellow_excitation_autofluorescence_percentile_95'), get_slope_wrapper(min_age,max_age,'green_yellow_excitation_autofluorescence_percentile_95')])
			except:	
				continue	
			color_vals = colorize.scale(lifespans, output_max=1)
			colors = colorize.color_map(color_vals, uint8=False)	
			plt.scatter(lifespans,results.y_est,c=colors)
			(m,b) = numpy.polyfit(lifespans, results.y_est, 1)
			yp=numpy.polyval([m,b], lifespans)
			mew, mewmew=zip(*sorted(zip(lifespans, yp)))
			mew, mewmew = (list(t) for t in zip(*sorted(zip(mew, mewmew))))
			plt.plot([mew[0],mew[-1]],[mewmew[0],mewmew[-1]],c='gray',alpha=.7)
			ftext="$r^{2}$ = "+(str)(round(results.R2,3))
			more_text="n= " + (str)(len(lifespans))
			plt.title('Regression on Slope and Mean ' +miRNA+'::GFP expression profiles '+(str)(min_age) + ' to '+(str)(max_age)+ ' '+key,fontdict=title_fontdict)
			plt.xlabel('Actual lifespan (days)',fontdict=label_fontdict)
			plt.ylabel("Predicted lifespan (days)",fontdict=label_fontdict)
			plt.figtext(.15,.8,ftext,fontsize=20,ha='left')
			plt.figtext(.8, .2, more_text, fontsize=20, ha='left')
			print('Saving '+save_name)	
			plt.savefig(save_name)
			plt.gcf().clf
			plt.close()


def straightened_worm_images(path_to_parameters):

	header,parameters=datafile.read_delimited(path_to_parameters+'/image_parameters.tsv')
	parameters=list(parameters)[0]
	maximum=parameters[header.index('maximum')]
	maximum=((int)(maximum))
	minimum=parameters[header.index('minimum')]
	minimum=((int)(minimum))
	gamma=parameters[header.index('gamma')]
	gamma=((int)(gamma))
	exp_root=parameters[header.index('exp_root')]
	position_root=parameters[header.index('position_root')]
	labels=position_root.split('/')
	annotations=load_data.read_annotations(exp_root)[labels[-1]]
	p=pathlib.Path(position_root)
	gfp=sorted(p.glob('* gfp.png'))
	bf=sorted(p.glob('* bf.png'))
	gfp_images=[freeimage.read(path) for path in gfp]
	bf_images=[freeimage.read(path) for path in bf]
	position_info,timepoint_info=annotations
	timepoints=timepoint_info.keys()
	width_tcks=[timepoint_info[key]['pose'][1] for key in timepoints]
	center_tcks=[timepoint_info[key]['pose'][0] for key in timepoints]
	ages=[timepoint_info[key]['age'] for key in timepoints]
	ages=[round(age/24,2)for age in ages]
	bf_splines=[]
	gfp_splines=[]


	for i in range(0, len(bf_images)):
		bf_spline=worm_spline.to_worm_frame(bf_images[i],center_tcks[i],width_tcks[i],width_margin=0,standard_length=1000)
		gfp_spline=worm_spline.to_worm_frame(gfp_images[i],center_tcks[i],width_tcks[i],width_margin=0,standard_length=1000)
		bf_splines.append(bf_spline)
		gfp_splines.append(gfp_spline)
	for bf_spline,gfp_spline,width_tck,age in zip(bf_splines,gfp_splines,width_tcks,ages):
		save_dir='/Volumes/9karray/Kinser_Holly/all_figures_figure1c/'+labels[4]+'/'+labels[5]
		if pathlib.Path(save_dir).exists()==False:
			os.makedirs(save_dir)
		mask=worm_spline.worm_frame_mask(width_tck,gfp_spline.shape,antialias=True)
		
		gfp_scaled=colorize.scale(gfp_spline,minimum,maximum,gamma=gamma).astype('uint8')
		bf_scaled=colorize.scale(bf_spline,1000,20000,gamma=1).astype('uint8')
	
		
		freeimage.write(numpy.dstack((bf_scaled,bf_scaled,bf_scaled,mask)),save_dir+'/'+(str)(age) +'_bf.png')
		freeimage.write(numpy.dstack((gfp_scaled,gfp_scaled,gfp_scaled,mask)),save_dir+'/'+(str)(age) +'_gfp.png')

def get_average_wrapper(min_age,max_age,feature):
	def get_average(worm):
		data=worm.get_time_range(feature+'_z',min_age,max_age)

		return numpy.mean(data[1])
	return get_average	
def get_slope_wrapper(min_age,max_age,feature):
	def get_slope(worm):
		data=worm.get_time_range(feature+'_z',min_age,max_age)
		slope, intercept, r_value, p_value, std_err=scipy.stats.linregress(data)

		return slope
	return get_slope	
#z transform all files in a directory and save out file with the same name (will keep old measurements)
def z_scoring(root_folder, age_feature='age',gfp_measures=gfp_measures,overwrite=True):
	subdirectories=glob(root_folder+'/*')
	for subdirectory in subdirectories:

		exp_dirs=glob(subdirectory+'/*')
		for exp_dir in exp_dirs:
			worms=process_worms([exp_dir], prefixes=[''])
			save_dir=pathlib.Path(exp_dir)
			if overwrite==False:
				continue
			for key,value in gfp_measures.items():
	
				try:
					worms.z_transform(key,age_feature=age_feature)
				except:
					continue	
			print('done transforming '+exp_dir)
			worms.write_timecourse_data(save_dir)			

			
def process_worms(paths,prefixes=''):
	wormies=[]
	for path, prefix in zip(paths,prefixes):
		worms=worm_data.read_worms(path+'/*.tsv',name_prefix=prefix)
		wormies=wormies+worms
	
	return wormies