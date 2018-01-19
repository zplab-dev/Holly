"""

@author: Holly
"""


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42



import numpy
import scipy.stats
import time
import pandas

import wormPhysiology.analyzeHealth.selectData as sd
import wormPhysiology.basicOperations.folderStuff as fs
import wormPhysiology.analyzeHealth.characterizeTrajectories as ct
import wormPhysiology.analyzeHealth.computeStatistics as cs

rainbow_colors=['purple','red','darkorange','gold']
other_colors=['indianred','darkred','hotpink','gold','indigo']


terms={
		'expressionarea_gfp': 'expression area', 
		'expressionareafraction_gfp': 'expression area fraction', 
		'expressionmean_gfp': 'expression area mean', 
		'highexpressionarea_gfp': 'high expression area', 
		'highexpressionareafraction_gfp': 'high expression area fraction', 
		'highexpressionintegrated_gfp': 'high expression area summed', 
		'highexpressionmean_gfp': 'high expression area mean', 
		'integrated99_gfp':'summed 99th percentile intensity',
		'integrated_gfp': 'summed intensity', 
		'intensity_70': '70th percentile intensity', 
		'intensity_80': '80th percentile intensity', 
		'intensity_90': '90th percentile intensity', 
		'intensity_95': '95th percentile intensity',
		'max_gfp': 'maximum intensity',
		'median_gfp': 'median intensity', 
		'percentile95_gfp': '95th percentile intensity', 
		'percentile99_gfp': '99th percentile intensity',
		'percentile99area_gfp':'99th percentile area',
		'total_size': 'total size',
		'health': 'health',
		'movement': 'movement',
		'autofluorescence': 'autofluorescence',
		'size': 'size'}

def control_for(measure, control_for, complete_worm_df, first_day, last_day, rescale=False):
	control_means, days, lifespans, worms, deleted_worms=get_data(complete_worm_df, [measure], first_day, last_day,rescale)

	means, days, lifespans, worms, deleted_worms=get_data(complete_worm_df,measure,first_day,last_day)
	slope_x, intercept_x,rvalue_x,pvalue_x,stderror_x=scipy.stats.linregress(x=control_means,y=means)
	slope_y,intercept_y,rvalue_y,pvalue_y,stderror_y=scipy.stats.linregress(x=control_means,y=lifespans)
	obsvalues_x=means
	predvalues_x=control*slope_x+intercept_x
	residuals_x=predvalues_x-obsvalues_x
	obsvalues_y=lifespans
	predvalues_y=control*slope_y+intercept_y
	residuals_y=predvalues_y-obsvalues_y

	pearson,spearman,yp=run_stats(residuals_x, residuals_y)
	
	return residuals_x, residuals_y, yp, correlation, spearman_correlation

def make_aggregate_plot(measures,first_day,last_day,miRNA,complete_worm_df, num_measures=4,save=True):

	all_correlations, days, all_p_values=get_all_correlations(measures, first_day, last_day, miRNA, complete_worm_df, save)
	ranked_correlations=[numpy.sum(all_correlations[i]) for i in range(0,len(all_correlations))]
	all_correlations,measures, ranked_correlations = (list(x) for x in zip(*sorted(zip(all_correlations,measures, ranked_correlations), reverse=True)))

	patches=[]
	for i in range(0, num_measures):
		patches.append(mlines.Line2D([],[],color=other_colors[i],label=terms[measures[i]],marker=',',linewidth=4.0))
	
	
	for i in range(0,num_measures):
		plt.plot(days, all_correlations[i], other_colors[i],linestyle='--',linewidth=2.0)
		plt.scatter(days,all_correlations[i],c=other_colors[i],marker='o',s=50,edgecolor=other_colors[i])
	
	plt.legend(handles=patches,loc='upper right', bbox_to_anchor=(0.35, 0.9),prop={'size':14,'family':'calibri'},frameon=True,fancybox=True,shadow=True,framealpha=.9,borderpad=1,labelspacing=.5,borderaxespad=1)
	plt.title('Correlation of '+miRNA +'::GFP expression with lifespan vs. timepoint measured', y=1.05,fontdict={'size':16, 'family':'calibri','weight':'bold'})
	ymin, ymax = plt.ylim()
	plt.ylim(ymin,ymax+.02)
	xmin,xmax=plt.xlim()
	plt.xlim(days[0]-.5,days[-1]+.5)
	plt.xlabel('Timepoint (day post-hatch)', fontdict={'size':14,'family':'calibri'},labelpad=10)
	plt.ylabel('Coefficient of determination ('+"$r^{2}$"+')', fontdict={'size':14,'family':'calibri'},labelpad=10)
	axis=plt.gca()
	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)
	fig_size = plt.rcParams["figure.figsize"]
	fig_size[0] = 15
	fig_size[1] = 9
	plt.rcParams["figure.figsize"] = fig_size
	plt.style.use('seaborn-white')
	plt.gcf()
	plt.savefig(complete_worm_df.save_directory+'/'+'top '+ (str)(num_measures)+' correlationtrajectories.pdf')	
	plt.show()
	plt.gcf().clf		
	

def get_all_correlations(measures,first_day,last_day,miRNA,complete_worm_df, save=True):
	
	all_correlations=[]
	all_p_values=[]
	for i in range(0, len(measures)):
		correlations,days,p_values=(get_correlations([measures[i]],first_day,last_day,miRNA,complete_worm_df))
		all_correlations.append(correlations)
		all_p_values.append(p_values)
	if save:
		dataframe=pandas.DataFrame(data=all_correlations,index=measures,columns=days)
		dataframe.to_csv(complete_worm_df.save_directory+'correlations.csv')

	return all_correlations,days,all_p_values		

def get_correlations(measure,first_day,last_day,miRNA, complete_worm_df):
	means, days, lifespans, worms, deleted_worms = get_data(complete_worm_df,measure,first_day,last_day)

	correlations=[]
	p_values=[]
	
	for i in range(0,len(means)):
		pearson,spearman,yp=run_stats(means[i],lifespans)
		correlations.append(pearson[0]**2)
		p_values.append(pearson[1])

	return correlations,days,p_values

def plot_all_correlations(first_day,last_day,measures,miRNA,complete_worm_df,show=False):

	for i in range(0,len(measures)):
		plot_correlations(first_day,last_day,[measures[i]],miRNA,complete_worm_df,show)


def plot_correlations(first_day,last_day,measure, miRNA, complete_worm_df,show=False):
	correlations,days,p_values =get_correlations(measure,first_day,last_day,miRNA,complete_worm_df)
	plt.style.use('seaborn-white')
	plt.scatter(days,correlations,c='indigo',marker='o',s=50,edgecolor='indigo')
	plt.plot(days,correlations, c='indigo',linewidth=2.0,linestyle='--')
	
	p_x=[]
	p_y=[]
	for i in range(0,len(p_values)):
		if p_values[i]<.05:
			p_y.append(correlations[i]+.03)
			p_x.append(days[i])
	plt.scatter(p_x,p_y,marker=(6,2,0),color='indigo',s=50)		
	plt.title('Correlation of '+miRNA + '::GFP expression ' + '('+ terms[measure[0]]+')'+ ' with lifespan vs. timepoint', y=1.05, fontdict={'size':16,'family':'calibri','weight':'bold'})
	plt.xlabel('Timepoint (day post-hatch)', fontdict={'size':14,'family':'calibri'})
	plt.ylabel('Coefficient of determination ('+"$r^{2}$"+')',fontdict={'size':14,'family':'calibri'})
	axis=plt.gca()
	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)
	ymin, ymax = plt.ylim()
	plt.ylim(ymin,ymax+.02)
	xmin,xmax=plt.xlim()
	plt.xlim(days[0]-.5,days[-1]+.5)
	fig_size = plt.rcParams["figure.figsize"]
	fig_size[0] = 15
	fig_size[1] = 9
	plt.rcParams["figure.figsize"] = fig_size
	plt.gcf()
	plt.savefig(complete_worm_df.save_directory+'/'+measure[0]+'correlationtrajectories.pdf')	
	if show:
		plt.show(block=False)
		time.sleep(2)
		plt.close()
	plt.gcf().clf		

def make_scatter_plots(all_measures,first_day,last_day, miRNA, complete_worm_df,show=False,control_for=None):

	all_means=[]
	
	for i in range(0, len(all_measures)):
		means, days, lifespans, worms, deleted_worms=get_data(complete_worm_df, [all_measures[i]], first_day, last_day)
		all_means.append(means)

	for i in range(0, len(all_means)):
		for j in range(0, len(days)):
			make_scatter_plot(all_means[i][j],lifespans,all_measures[i],days[j],miRNA,complete_worm_df,show)


def make_scatter_plot(means, lifespans, measure, day, miRNA, complete_worm_df,show=False,control_for=None):

	life_cohorts, bin_lifes, my_bins, my_colors=sd.life_cohort_bins(complete_worm_df)

	pearson,spearman,yp=run_stats(means,lifespans)
	
	colors=[]
	
	for i in range(0,len(lifespans)):
		for j in range(0,len(my_bins)):
			if lifespans[i]<my_bins[j][1] and lifespans[i]>my_bins[j][0]:
				colors.append(rainbow_colors[j])		

	plt.scatter(means,lifespans,c=colors, edgecolors=colors,s=35,alpha=.7)
	plt.plot(means,yp, c='gray')
	plt.style.use('seaborn-white')

	title=''
	y_label=''
	save_name=''
	if measure.endswith('gfp'):
		title="Mean "+miRNA+ "::GFP expression at "+ (str)(day) + " dph vs. lifespan"
		y_label="Mean "+miRNA + "::GFP expression ( "+terms[measure]+ " )"
		save_name= "Mean "+miRNA + " expression "+ (str)(day)+ " dph " + measure
	if measure.endswith('size'):
		title="Total size at "+ (str)(day) + " dph vs. lifespan"
		y_label="Total size (pixels)"
		save_name= "Total size  "+(str)(day)+ " dph " 
	if measure.startswith('intensity'):
		title="Mean autofluorescence at "+ (str)(day) + " dph vs. lifespan"
		y_label="Mean autofluorescence ( "+measures[i]+ " )"
		save_name= "Mean autofluorescence "+(str)(day)+ " dph " + measure 		

	
	plt.title(title,y=1.05,fontdict={'size':16,'family':'calibri','weight':'bold'})
	plt.xlabel("Lifespan (days)",fontdict={'size':14,'family':'calibri'})
	plt.ylabel(y_label,fontdict={'size':14,'family':'calibri'})
	
	axis=plt.gca()
	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)
	
	if pearson[1]<.00001:
		p="p<.00001"
	else:
		p="p=" + ''+ (str)(round(pearson[1],3))
	if spearman[1]<.00001:
		spearman_p="p<.00001"
	else:
		spearman_p="p=" + '' + (str)(round(spearman[1],3))
				
	ftext="$r^{2}$ ="+(str)(round(pearson[0]**2,3))+" "+p
	gtext="rho="+(str)(round(spearman[0],3))+" "+spearman_p
	more_text="n= " + (str)(len(means))
	
	plt.figtext(.5,.85,ftext,fontsize=14,ha='left', family='calibri')
	plt.figtext(.15,.8,gtext,fontsize=14,ha='left',family='calibri')
	plt.figtext(.8, .2, more_text, fontsize=14, ha='left',family='calibri')
	
	fig_size = plt.rcParams["figure.figsize"]
	fig_size[0] = 15
	fig_size[1] = 9
	plt.rcParams["figure.figsize"] = fig_size
	
	plt.gcf()
	plt.savefig(complete_worm_df.save_directory+'/'+save_name+'.pdf')
	if show:
		plt.show(block=False)
		time.sleep(1)
		plt.close()
	plt.gcf().clf


def run_stats(x_list,y_list):
	
	pearson=numpy.asarray(scipy.stats.pearsonr(x_list, y_list))
	spearman=numpy.asarray(scipy.stats.spearmanr(x_list, y_list))
	(m,b) = numpy.polyfit(x_list, y_list, 1)
	yp=numpy.polyval([m,b], x_list)
	
	return (pearson,spearman, yp)

def get_data(complete_worm_df,measure,first_day,last_day,rescale=False):

	select_times=complete_worm_df.times[complete_worm_df.times.index(first_day):(complete_worm_df.times.index(last_day)+8)]
	lifespans=sd.get_lifespans(complete_worm_df)/24
	worms=complete_worm_df.worms
	data=complete_worm_df.mloc(worms,measure,select_times)
	
	last_day=(int)((float)(last_day))
	first_day=(int)((float)(first_day))
	
	days=[]
	
	for k in range(first_day, last_day+1,1):
		days.append(k)


	datakins=[]
	
	for i in range(0, (last_day-first_day)+1):
		datakins.append((data[:,0,i*8:((i*8)+8)]))
	
	means=[]
	
	for i in range(0,len(datakins)):
		means.append([])
	for i in range(0, len(datakins)):
		for j in range(0,len(datakins[i])):
			means[i].append(numpy.mean(datakins[i][j]))

	nan_positions=[]        
	
	for j in range(0,len(means[0])):
		if numpy.isnan(means[len(means)-1][j]):
				nan_positions.append(j)
	for i in range(0, len(means)):	
		mew=numpy.delete(means[i],nan_positions)
		means[i]=mew
	
	if rescale:
		means=numpy.asarray(means)
		rescaled_means=numpy.empty(means.shape)
		for i in range(0,len(means)):
			median=numpy.median(means[i])
			for j in range(0,len(means[i])):
				rescaled_means[i][j]=means[i][j]/median
		means=rescaled_means

	lifespans=numpy.delete(lifespans,nan_positions)
	deleted_worms=[worms[position] for position in nan_positions]
	worms=numpy.delete(worms,nan_positions)
	
	return means, days, lifespans, worms, deleted_worms

def make_dataframe(data_directories, save_directory):
	work_dir='/Volumes/purplearray/Kinser_Holly/work_dir' 
	human_dir='/Volumes/purplearray/Kinser_Holly/utilities'
	extra_dirs=[]
	experiment_dirs=[]
	annotations_dirs=[]
	for i in range(0, len(data_directories)):
		extra_dirs.append(None)
		experiment_dirs.append(None)
		annotations_dirs.append(None)
	directory_bolus=fs.DirectoryBolus(work_dir, human_dir, data_directories,extra_dirs, experiment_dirs, annotations_dirs, done=len(data_directories), ready=len(data_directories))
	complete_worm_df=ct.CompleteWormDF(directory_bolus,save_directory) 
	return complete_worm_df


def measure_slopes(measure, complete_worm_df, first_day, last_day,miRNA):
	worms=complete_worm_df.worms
	lifespans=sd.get_lifespans(complete_worm_df)/24
	select_times=complete_worm_df.times[complete_worm_df.times.index(first_day):(complete_worm_df.times.index(last_day)+1)]
	data=complete_worm_df.mloc(worms,[measure],select_times)
	times=[float(i) for i in select_times]

	life_cohorts, bin_lifes, my_bins, my_colors=sd.life_cohort_bins(complete_worm_df)
	slopes=[]

	for i in range(0, len(data)):
		slope,intercept,r,p,std_err=scipy.stats.linregress(times, data[i][0])
		slopes.append(slope)
	slopes=numpy.asarray(slopes)
	lifespans=lifespans[numpy.isnan(slopes)==False]
	slopes=slopes[numpy.isnan(slopes)==False]
	
	colors=[]
	for i in range(0,len(slopes)):
		for j in range(0,len(life_cohorts)):
			if i in life_cohorts[j]:
				colors.append(rainbow_colors[j])

	plt.style.use('seaborn-white')
	axis=plt.gca()
	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)
	plt.title(miRNA +" slope (days "+(str)(first_day)+" to "+ (str)(last_day)+")" +" vs. lifespan",fontdict={'size':16,'family':'calibri','weight':'bold'})
	plt.xlabel("Slope of "+miRNA+ "::GFP ("+terms[measure]+")",fontdict={'size':14,'family':'calibri'})
	plt.ylabel("Lifespan (days)",fontdict={'size':14,'family':'calibri'})
	plt.scatter(slopes,lifespans,c=colors, edgecolors=colors,s=35,alpha=.7)
	correlation=scipy.stats.pearsonr(slopes, lifespans)
	correlation = numpy.asarray(correlation)
	(m,b) = numpy.polyfit(slopes, lifespans, 1)
	yp=numpy.polyval([m,b], slopes)
	plt.plot(slopes,yp,c='gray')
	if correlation[1]<.00001:
		p="p<.00001"
	else:
		p="p=" + ''+ (str)(round(correlation[1],3))
	ftext="$r^{2}$ ="+(str)(round(correlation[0]**2,3))+" "+p
	plt.figtext(.65,.8,ftext,fontsize=14,ha='left')
	
	fig_size = plt.rcParams["figure.figsize"]
	fig_size[0] = 15
	fig_size[1] = 9
	plt.rcParams["figure.figsize"] = fig_size

	plt.savefig(complete_worm_df.save_directory+'/'+'slope ' +(str)(first_day) +' to '+(str)(last_day)+' dph.pdf')	
	plt.show()
	plt.gcf().clf


def plot_survival_curve(complete_worm_df, save_directory):
	
	lifespans=sd.get_lifespans(complete_worm_df)/24

	max_life=lifespans.max()
	min_life=lifespans.min()

	days=numpy.arange(min_life-1,max_life+1,.5)
	
	percent_alive = []
	
	for i in days:
		count =0
		for item in lifespans:
			if item > i:
				count=count+1
		percent_alive.append((count/len(lifespans))*100)		
	
	plt.style.use('seaborn-white')
	plt.plot(days,percent_alive, color='goldenrod',linewidth=2.0)
	
	axis=plt.gca()
	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)
	plt.title("Survival curve",y=1.05,fontdict={'size':16,'family':'calibri','weight':'bold'})
	plt.xlabel("Lifespan (days)",fontdict={'size':12,'family':'calibri'})
	plt.ylabel("Survival (%)",fontdict={'size':12,'family':'calibri'})
	ftext="mean lifespan="+(str)(round(lifespans.mean(),1)) + " days"
	even_more_text="median lifespan = " + (str)(round(numpy.percentile(lifespans,50),1)) + " days"
	more_text="n= " + (str)(lifespans.size)
	plt.figtext(.5,.8,ftext,fontsize=12,ha='left',family='calibri')
	plt.figtext(.8, .2, more_text, fontsize=12, ha='left', family='calibri')
	plt.figtext(.5, .75, even_more_text, fontsize=12, ha='left',family='calibri')
	fig_size = plt.rcParams["figure.figsize"]
	fig_size[0] = 15
	fig_size[1] = 9
	plt.rcParams["figure.figsize"] = fig_size
	plt.gcf()
	plt.savefig(save_directory+'/survival_curve.pdf')
	plt.show()
	plt.gcf().clf

def plot_smoothed_trajectories(complete_worm_df, measure, save_directory,bin_width_days=2, miRNA=''):
	
	smoothed_data=complete_worm_df.smooth_trajectory(measure,extra_arguments={'smooth_iterations': 3, 'polyorder': 1, 'windowLength': 9, 'smooth_mode': 'savitzky_golay'})
	ages=complete_worm_df.ages
	life_cohorts,bin_lifes,my_bins,my_colors=sd.life_cohort_bins(complete_worm_df,bin_width_days=bin_width_days)
	averages=[]
	
	plt.style.use('seaborn-white')
	plt.xlabel('Age (days)', fontdict={'size':12,'family':'calibri'})
	axis=plt.gca()
	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)
	plt.title('Average '+miRNA+'::GFP profiles of different lifespan cohorts',y=1.05,fontdict={'size':14,'family':'calibri', 'weight':'bold'})
	plt.ylabel(miRNA+ '::GFP expression ('+ terms[measure]+')', fontdict={'size':12,'family':'calibri'})
	more_text='n =  '
	
	for i in range(0, len(life_cohorts)):
		more_text=more_text+'       '+ (str)(len(life_cohorts[i]))

	plt.figtext(.5, .2, more_text, fontsize=12, ha='left',family='calibri')
	
	for i in range(0, len(life_cohorts)):
		cohort_data = smoothed_data[life_cohorts[i]]
		cohort_data = cohort_data[~numpy.isnan(cohort_data).all(axis = 1)]
		data=numpy.mean(cohort_data, axis = 0)
		mew=numpy.copy(data)
		mew[ages.index(round(bin_lifes[i]),1):]='NaN'
		averages.append(mew)
		
	for i in range(0, len(life_cohorts)):
		plt.plot(ages, averages[i], color=rainbow_colors[i],linewidth=2.0)
	save_name='average'+measure+'profiles'	

	fig_size = plt.rcParams["figure.figsize"]
	fig_size[0] = 15
	fig_size[1] = 9
	plt.rcParams["figure.figsize"] = fig_size

	plt.gcf()
	plt.savefig(save_directory+'/'+save_name+'.pdf')
	plt.show()

def plot_trajectories(complete_worm_df, measure, save_directory,bin_width_days=2, miRNA=''):
	ages=complete_worm_df.ages
	life_cohorts,bin_lifes,my_bins,my_colors=sd.life_cohort_bins(complete_worm_df,bin_width_days=bin_width_days)
	averages=[]
	
	plt.style.use('seaborn-white')
	plt.xlabel('Age (days)', fontdict={'size':12,'family':'calibri'})
	axis=plt.gca()
	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)
	plt.title('Average '+miRNA+'::GFP profiles of different lifespan cohorts',y=1.05,fontdict={'size':14,'family':'calibri', 'weight':'bold'})
	plt.ylabel(miRNA+ '::GFP expression ('+ terms[measure]+')', fontdict={'size':12,'family':'calibri'})
	more_text='n =  '
	
	for i in range(0, len(life_cohorts)):
		more_text=more_text+'       '+ (str)(len(life_cohorts[i]))

	plt.figtext(.5, .2, more_text, fontsize=12, ha='left',family='calibri')
	
	for i in range(0, len(life_cohorts)):
		averages.append(sd.cohort_trace(complete_worm_df,life_cohorts[i],measure,False))
	for i in range(0, len(life_cohorts)):
		plt.plot(ages, averages[i], color=rainbow_colors[i])
		
	save_name='averages'+measure+'profiles'	

	fig_size = plt.rcParams["figure.figsize"]
	fig_size[0] = 15
	fig_size[1] = 9
	plt.rcParams["figure.figsize"] = fig_size
	
	plt.gcf()
	plt.savefig(out_dir+'/'+save_name+'.pdf')
	plt.show()
	plt.gcf().clf		

def make_cv_plots(cv_list,days,measure,miRNA,save_dir):
	plt.style.use('seaborn-white')
	plt.plot(days,cv_list, c='mediumseagreen')
	plt.title('Coefficient of variation '+miRNA + ' expression ' + '('+ terms[measure]+')'+ ' vs. timepoint', y=1.05, fontdict={'size':14,'family':'calibri'})
	plt.xlabel('Timepoint (dph)', fontdict={'size':12,'family':'calibri'})
	plt.ylabel('Coefficient of variation (mean/std)',fontdict={'size':12,'family':'calibri'})
	axis=plt.gca()
	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)
	ymin, ymax = plt.ylim()
	plt.ylim(ymin,ymax+.02)
	plt.gcf() 
	plt.savefig(save_dir+'/'+measure+'correlationtrajectories.png')	
	plt.show(block=False)
	time.sleep(1)
	plt.close()
	plt.gcf().clf


def plot_spread(means,measure,save_dir,miRNA):
	plt.style.use('seaborn-white')
	plt.hist(means, facecolor='deepskyblue', alpha=.5)
	plt.title('Histogram of '+miRNA + ' expression ' + '('+ terms[measure]+')', y=1.05, fontdict={'size':14,'family':'calibri'})
	plt.xlabel(miRNA + ' expression ' + '('+ terms[measure]+')', fontdict={'size':12,'family':'calibri'})
	plt.ylabel('Frequency',fontdict={'size':12,'family':'calibri'})
	axis=plt.gca()
	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)
	ymin, ymax = plt.ylim()
	xmin,xmax=plt.xlim()
	cv=numpy.std(means)/numpy.mean(means)
	plt.ylim(ymin,ymax+5)
	plt.xlim(xmin-numpy.std(means),xmax+numpy.std(means))
	more_text="cv= " + (str)(cv)
	plt.figtext(.2, .8, more_text, fontsize=13, ha='left',family='calibri')
	plt.gcf()

	plt.savefig(save_dir+'/'+measure+'histogram.png')	
	plt.show(block=False)
	time.sleep(1)
	plt.close()
	plt.gcf().clf
	return cv


