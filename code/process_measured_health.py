"""

@author: Holly
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
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










import numpy
import scipy.stats
import time
import pandas
import math

import wormPhysiology.analyzeHealth.selectData as sd
import wormPhysiology.basicOperations.folderStuff as fs
import wormPhysiology.analyzeHealth.characterizeTrajectories as ct
import wormPhysiology.analyzeHealth.computeStatistics as cs
rainbow_colors=['indigo','darkorchid','mediumvioletred','crimson','red','orangered','darkorange','orange','gold','yellow']
#rainbow_colors=['magenta','purple','red','darkorange','gold','yellow','cyan','pink','black','grey']
other_colors=['indianred','darkred','hotpink','gold','indigo']

color_terms={'percentile95_gfp': 'indigo',
			'integrated_gfp': 'darkred',
			'median_gfp': 'hotpink',
			'max_gfp': 'gold',
			'percentile99_gfp': 'indianred',
			'expressionarea_gfp': 'darkorchid',
			'highexpressionintegrated_gfp': 'mediumvioletred',
			'expressionareafraction_gfp':'darkturquoise',
			'highexpressionarea_gfp': 'grey',
			'expressionarea_gfp': 'crimson',
			'percentile99area_gfp': 'black',
			'highexpressionareafraction_gfp':'red',
			'integrated99_gfp': 'mediumpurple',
			'total_size':'forestgreen',
			'expressionmean_gfp': 'teal',
			'highexpressionmean_gfp':'brown'	
	
}
terms={
		'expressionarea_gfp': 'expression area',
		'expressionarea_rfp': 'expression area', 
		'expressionareafraction_gfp': 'expression area fraction',
		'expressionareafraction_rfp': 'expression area fraction', 
		'expressionmean_gfp': 'expression area mean',
		'expressionmean_rfp': 'expression mean', 
		'highexpressionarea_gfp': 'high expression area', 
		'highexpressionarea_rfp': 'high expression area',
		'highexpressionareafraction_gfp': 'high expression area fraction', 
		'highexpressionareafraction_rfp': 'high expression area fraction',
		'highexpressionintegrated_gfp': 'high expression area summed',
		'highexpressionintegrated_rfp': 'high expression area summed', 
		'highexpressionmean_gfp': 'high expression area mean', 
		'highexpressionmean_rfp': 'high expression area mean',
		'integrated99_gfp':'summed 99th percentile intensity',
		'integrated99_rfp': 'summed 99th percentile intensity',
		'integrated_gfp': 'summed intensity',
		'integrated_rfp': 'summed intensity', 
		'intensity_70': 'autofluorescence (70th percentile intensity)', 
		'intensity_80': 'autofluorescence (80th percentile intensity)', 
		'intensity_90': 'autofluorescence (90th percentile intensity)', 
		'intensity_95': 'autofluorescence (95th percentile intensity)',
		'max_gfp': 'maximum intensity',
		'max_rfp': 'maximum intensity',
		'median_gfp': 'median intensity', 
		'median_rfp': 'median intensity',
		'percentile95_gfp': '95th percentile intensity',
		'percentile95_rfp': '95th percentile intensity', 
		'percentile99_gfp': '99th percentile intensity',
		'percentile99_rfp': '99th percentile intensity',
		'percentile99area_gfp':'99th percentile area',
		'percentile99area_rfp': '99th percentile area',
		'total_size': 'total size',
		'health': 'health',
		'movement': 'movement',
		'autofluorescence': 'autofluorescence',
		'total_size': 'size'}

def cohort_survival_curve(df, cohort_info=None,make_labels=True,bin_width_days=1):
	
	lifespans = sd.get_lifespans(df)/24
	life_histogram = numpy.histogram(lifespans, density = True, bins = 1000)
	life_times = life_histogram[1]
	cumulative_death = life_histogram[0]/numpy.sum(life_histogram[0])
	cumulative_death = numpy.cumsum(cumulative_death)
	cumulative_life = 1 - cumulative_death


	cumulative_life = numpy.insert(cumulative_life, 0, 1)
	cumulative_life = numpy.insert(cumulative_life, 0, 1)
	life_times = numpy.insert(life_times, 0, 0)
	
	if cohort_info is None:
		(life_cohorts, bin_lifes, my_bins, my_colors)=sd.life_cohort_bins(df,bin_width_days=bin_width_days)
	else:
		(life_cohorts, bin_lifes, my_bins, my_colors) = cohort_info
	

	cohort_lifes = numpy.array([lifespans[a_cohort] for a_cohort in life_cohorts])
	cohort_mins = [numpy.min(cohort_life) for cohort_life in cohort_lifes]
	cohort_mins = [numpy.argmin(numpy.abs(life_times - cohort_min)) for cohort_min in cohort_mins]  

	
	plt.plot(life_times[:cohort_mins[0] + 1], cumulative_life[:cohort_mins[0] + 1], color = 'black')
	for i in range(0, len(cohort_mins) - 1):
			
		plt.plot(life_times[cohort_mins[i]: cohort_mins[i+1] + 1], cumulative_life[cohort_mins[i]: cohort_mins[i+1] + 1], color = my_colors[i])
		
		plt.plot(life_times[cohort_mins[-1]:], cumulative_life[cohort_mins[-1]:], color = my_colors[-1])
	plt.style.use('seaborn-white')
	
	plt.title('Survival curve', y=1.05, fontdict={'size':26,'weight':'bold', 'family':'calibri'})
	plt.xlabel('Lifespan (days)', fontdict={'size':20,'family':'calibri'})
	plt.ylabel('Survival (%)',fontdict={'size':20,'family':'calibri'})
	plt.xlim([3,lifespans.max()//1+1])
	plt.ylim([0.0, 1.1])

	ftext="mean lifespan="+(str)(round(lifespans.mean(),1)) + " days"
	even_more_text="median lifespan = " + (str)(round(numpy.percentile(lifespans,50),1)) + " days"
	more_text="n= " + (str)(lifespans.size)
	plt.figtext(.7,.8,ftext,fontsize=20,ha='left')
	plt.figtext(.8, .2, more_text, fontsize=20, ha='left')
	plt.figtext(.7, .75, even_more_text, fontsize=20, ha='left')


	plt.gcf()
	plt.savefig(df.save_directory+'/Survival_curve.pdf')
	plt.show()
	plt.gcf().clf
def cohort_histogram(df,cohort_info=None,make_labels=True,bin_width_days=2):
	lifespans = sd.get_lifespans(df)/24
	life_histogram = numpy.histogram(lifespans, density = True, bins = 1000)
	life_times = life_histogram[1]
	cumulative_death = life_histogram[0]/numpy.sum(life_histogram[0])
	cumulative_death = numpy.cumsum(cumulative_death)
	cumulative_life = 1 - cumulative_death

	# Insert a point to start with zero days post hatch.
	cumulative_life = numpy.insert(cumulative_life, 0, 1)
	cumulative_life = numpy.insert(cumulative_life, 0, 1)
	life_times = numpy.insert(life_times, 0, 0)

	# Figure out where to start each line.
	if cohort_info is None:
		(life_cohorts, bin_lifes, my_bins, my_colors)=sd.life_cohort_bins(df,bin_width_days=bin_width_days)
	else:
		(life_cohorts, bin_lifes, my_bins, my_colors) = cohort_info
	

	cohort_lifes = numpy.array([lifespans[a_cohort] for a_cohort in life_cohorts])
	cohort_mins = [numpy.min(cohort_life) for cohort_life in cohort_lifes]
	cohort_mins = [numpy.argmin(numpy.abs(life_times - cohort_min)) for cohort_min in cohort_mins]  

	bin_width=2
	plt.style.use('seaborn-white')
	
	plt.title('Cohorts across distribution of lifespans', y=1.05,fontdict={'size':26, 'family':'calibri', 'weight':'bold'})
	plt.xlabel('Lifespan (days)', fontdict={'size':20,'family':'calibri'})
	plt.ylabel('Number of individuals',fontdict={'size':20,'family':'calibri'})
	plt.xlim(lifespans.min()-2,lifespans.max()+2)
	
	#more_text='n =  '
	#for i in range(0, len(life_cohorts)):
		#more_text=more_text+'       '+ (str)(len(life_cohorts[i]))

	#plt.figtext(.5, .2, more_text, fontsize=20, ha='center',family='calibri',weight='bold')
	
	for i in range(0, len(cohort_lifes)):
		my_range = (numpy.floor(numpy.min(cohort_lifes[i])), numpy.ceil(numpy.max(cohort_lifes[i])))
		(n, bins, patches) = plt.hist(cohort_lifes[i], range = my_range, bins = 2/bin_width, normed = False, facecolor = rainbow_colors[i], alpha = .7, linewidth = 0)
	


	plt.savefig(df.save_directory+'/'+'lifespanhist.pdf')	
	plt.show()
	plt.gcf().clf	

def consistency_index (measure, complete_worm_df, first_day, last_day):
	means, days, lifespans, worms, deleted_worms = get_data(complete_worm_df,[measure],first_day,last_day)
	medians=numpy.median(means,axis=1)
	consistencies=numpy.zeros(worms.shape)
	for i in range(0, len(consistencies)):
		num_higher=1
		num_lower=1
		for j in range(0,len(days)):
			if means[j][i]>medians[j]:
				num_higher+=1
			else:
				num_lower+=1
		index=math.log(num_higher/num_lower,2)
		consistencies[i]=index
	
	null_means=[]
	
	for i in range(0, len(days)):
		null_means.append(numpy.random.choice(means[i],1000))

	null_consistencies=numpy.zeros(1000)
	null_medians=numpy.median(null_means,axis=1)
	for i in range(0, len(null_consistencies)):
		num_higher=1
		num_lower=1
		for j in range(0, len(days)):
			if null_means[j][i]>null_medians[j]:
				num_higher+=1
			else:
				num_lower+=1
		index=math.log(num_higher/num_lower,2)
		null_consistencies[i]=index

	p_values=[]

	for i in range(0,len(consistencies)):
		count=0
		for j in range(0, len(null_consistencies)):
			if null_consistencies[j]>=consistencies[i] and consistencies[i]>0:
				count+=1
			if null_consistencies[j]<=consistencies[i] and consistencies[i]<0:
				count+=1
			if consistencies[i]==0:
				count+=100
		p_values.append(count/1000)		
	
	colors=[]
	for i in range(0, len(consistencies)):
		if p_values[i]<.05:
			colors.append('red')
		else: 
			colors.append('grey')		
	x_positions=numpy.random.sample(consistencies.shape)
	x_positions_null=numpy.random.sample(null_consistencies.shape)+3
	
	plt.style.use('seaborn-white')
	
	plt.scatter(x_positions,consistencies, color=colors)
	plt.scatter(x_positions_null,null_consistencies,color='black')
	plt.ylim(-5,5)
	ax = plt.gca()
	ax.set_xticks([.5,3.5])
	ax.set_xticklabels(['Real','Shuffled'],fontdict={'size':14,'family':'calibri'})
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.axhline(y=0,color='black',linestyle='dashed')

	plt.show()	
		
			
	return consistencies, p_values				
def correlate_two_measures(measure1, measure2, first_day, last_day,miRNA,complete_worm_df):

	means1, days, lifespans, worms, deleted_worms = get_data(complete_worm_df,[measure1],first_day,last_day)
	means2, days, lifespans, worms, deleted_worms = get_data(complete_worm_df, [measure2],first_day,last_day)

	life_cohorts, bin_lifes, my_bins, my_colors=sd.life_cohort_bins(complete_worm_df)

	colors=[]
	for i in range(0,len(lifespans)):
		for j in range(0,len(my_bins)):
			if lifespans[i]<my_bins[j][1] and lifespans[i]>my_bins[j][0]:
				colors.append(rainbow_colors[j])		
	
	plt.style.use('seaborn-white')
	
	for i in range(0, len(days)):
		pearson,spearman,yp=run_stats(means1[i],means2[i])				
		plt.scatter(means1[i],means2[i],c=colors, edgecolors=colors,s=35,alpha=.7)
		plt.plot(means1[i],yp, c='gray')
		title="Mean "+miRNA+ "::GFP expression vs mean "+ terms[measure2]+" at " +(str)(days[i]) + " dph"
		x_label="Mean "+miRNA + "::GFP expression ( "+terms[measure1]+ " )"
		y_label="Mean "+ terms[measure2]
		save_name= "Mean "+miRNA + " expression "+ (str)(days[i])+ " dph " + measure1+measure2
		plt.title(title,y=1.05,fontdict={'size':26,'weight':'bold','family':'calibri'})
		plt.ylabel(y_label,fontdict={'size':20,'family':'calibri'})
		plt.xlabel(x_label,fontdict={'size':20,'family':'calibri'})
	
	
		if pearson[1]<.00001:
			p="p<.00001"
		else:
			p="p=" + ''+ (str)(round(pearson[1],3))
		if spearman[1]<.00001:
			spearman_p="p<.00001"
		else:
			spearman_p="p=" + '' + (str)(round(spearman[1],3))
				
		ftext="$r^{2}$ =" +(str)(round(pearson[0]**2,3))+" "+p
		gtext=r'$\rho$'+" = "+(str)(round(spearman[0],3))+" "+spearman_p
		more_text="n= " + (str)(len(means1[i]))
	
		plt.figtext(.15,.85,ftext,fontsize=20,ha='left')
		plt.figtext(.15,.8,gtext,fontsize=20,ha='left')
		plt.figtext(.8, .2, more_text, fontsize=20, ha='left')
	
		
	
		plt.gcf()
		plt.savefig(complete_worm_df.save_directory+'/'+save_name+'.pdf')
		plt.show(block=False)
		time.sleep(3)
		plt.close()
		plt.gcf().clf

def plot_all_CVs(first_day,last_day,measures,miRNA,complete_worm_df,show=False, rescale=[]):

	for i in range(0,len(measures)):
		plot_CVs(first_day,last_day,[measures[i]],miRNA,complete_worm_df,show,rescale)

def plot_CVs(first_day,last_day,measure, miRNA, complete_worm_df,show=False,rescale=[]):
	cvs,days=get_CVs(measure,first_day,last_day,miRNA,complete_worm_df,rescale)
	plt.style.use('seaborn-white')
	plt.scatter(days,cvs,c='indigo',marker='o',s=50,edgecolor='indigo')
	plt.plot(days,cvs, c='indigo',linestyle='--')
	plt.title('Coefficient of variation of '+miRNA + '::GFP expression ' + '('+ terms[measure[0]]+')'+ ' vs. timepoint', y=1.05, fontdict={'size':26,'weight':'bold','family':'calibri'})
	plt.xlabel('Timepoint (day post-hatch)', fontdict={'size':20,'family':'calibri'})
	plt.ylabel('Coefficient of variation',fontdict={'size':20,'family':'calibri'})

	ymin, ymax = plt.ylim()

	xmin,xmax=plt.xlim()
	plt.xlim(days[0]-.5,days[-1]+.5)

	plt.gcf()
	
	plt.savefig(complete_worm_df.save_directory+'/'+'CVtrajectories'+measure[0]+'.pdf')	
	if show:
		plt.show(block=False)
		time.sleep(2)
		plt.close()
	plt.gcf().clf			

def get_CVs(measure,first_day,last_day,miRNA, complete_worm_df,rescale=[]):

	if len(rescale)==0:
		means, days, lifespans, worms, deleted_worms = get_data(complete_worm_df,measure,first_day,last_day)
	else:
		means, days, lifespans, worms, deleted_worms = run_rescaling_multiple(rescale,measure,first_day,last_day)	
	
	cvs=[]

	for i in range(0,len(means)):
		cv=round(numpy.std(means[i])/numpy.mean(means[i]),2)
		cvs.append(cv)

	return cvs,days

def get_all_CVs(measures,first_day,last_day,miRNA,complete_worm_df, save=True,rescale=[]):
	
	all_cvs=[]
	
	for i in range(0, len(measures)):
		cvs,days=(get_CVs([measures[i]],first_day,last_day,miRNA,complete_worm_df,control_for,rescale))
		all_cvs.append(cvs)
		
	if save:
		dataframe=pandas.DataFrame(data=all_cvs,index=measures,columns=days)
		dataframe.to_csv(complete_worm_df.save_directory+'/cvs.csv')

	return all_correlations,days

def plot_all_spreads(all_measures,first_day,last_day, miRNA, complete_worm_df,show=False,rescale=[]):
	all_means=[]
	
	for i in range(0, len(all_measures)):
		if len(rescale)==0:
			means, days, lifespans, worms, deleted_worms=get_data(complete_worm_df, [all_measures[i]], first_day, last_day)
		if len(rescale)>0:
			means, days, lifespans, worms, deleted_worms=run_rescaling_multiple(rescale, [all_measures[i]], first_day, last_day)	
		
		all_means.append(means)
	
	for i in range(0, len(all_means)):
		for j in range(0, len(days)):
			plot_spread(all_means[i][j],all_measures[i],complete_worm_df,days[j],miRNA,show)
	
	print(deleted_worms)
	print("deleted "+(str)(round(100*(len(deleted_worms)/len(worms)),1))+" percent of worms")

def plot_spread(means,measure,complete_worm_df,day, miRNA,show=False):
	plt.style.use('seaborn-white')
	plt.hist(means, facecolor='indigo')
	plt.title('Distribution of '+miRNA + ' expression ' + '('+ terms[measure]+')'+' at '+ (str)(day) + ' dph', y=1.05, fontdict={'size':26,'weight':'bold','family':'calibri'})
	plt.xlabel(miRNA + ' expression ' + '('+ terms[measure]+')', fontdict={'size':20,'family':'calibri'})
	plt.ylabel('Frequency',fontdict={'size':20,'family':'calibri'})

	ymin, ymax = plt.ylim()
	xmin,xmax=plt.xlim()
	
	cv=round(numpy.std(means)/numpy.mean(means),2)
	plt.ylim(ymin,ymax+5)
	plt.xlim(xmin-numpy.std(means),xmax+numpy.std(means))
	more_text="coefficient of variation = " + (str)(cv)
	plt.figtext(.2, .8, more_text, fontsize=20, ha='left')
	plt.gcf()

	plt.savefig(complete_worm_df.save_directory+'/'+'Histogram '+measure+" at day "+(str)(day)+'.pdf')
	if show:
		plt.show(block=False)
		time.sleep(1)
		plt.close()
	plt.gcf().clf		
	
	

def make_aggregate_plot(measures,first_day,last_day,miRNA,complete_worm_df, num_measures=4,save=True,control_for='',rescale=[]):

	all_correlations, days, all_p_values=get_all_correlations(measures, first_day, last_day, miRNA, complete_worm_df, save,control_for,rescale)
	summed_correlations=[numpy.sum(all_correlations[i]) for i in range(0,len(all_correlations))]
	ranked_correlations=sorted(summed_correlations, reverse=True)

	ranked_measures=[measures[summed_correlations.index(ranked_correlations[i])] for i in range(0,num_measures)]
	extra_ranked_correlations=[all_correlations[summed_correlations.index(ranked_correlations[i])] for i in range(0,num_measures)]
	
	patches=[]
	for i in range(0, len(ranked_measures)):
		patches.append(mlines.Line2D([],[],color=color_terms[ranked_measures[i]],label=terms[ranked_measures[i]],marker=','))
	
	
	for i in range(0,len(ranked_measures)):
		plt.plot(days, extra_ranked_correlations[i], color=color_terms[ranked_measures[i]],linestyle='--')
		plt.scatter(days,extra_ranked_correlations[i],c=color_terms[ranked_measures[i]],marker='o',s=50,edgecolor=color_terms[ranked_measures[i]])
	
	plt.legend(handles=patches,loc='best', prop={'size':18},labelspacing=.5)
	
	if len(control_for)>0:
		plt.title('Correlation of '+miRNA +' expression with lifespan vs. timepoint measured, controlled for '+terms[control_for], y=1.05,fontdict={'size':26, 'weight':'bold','family':'calibri'})
	else:
		plt.title('Correlation of '+miRNA +' expression with lifespan vs. timepoint measured', y=1.05,fontdict={'size':26,'weight':'bold','family':'calibri'})
	
	ymin, ymax = plt.ylim()
	plt.ylim(0,.5)
	#plt.ylim(ymin,ymax+.1)
	
	xmin,xmax=plt.xlim()
	plt.xlim(days[0]-.5,days[-1]+.5)

	plt.xlabel('Timepoint (day post-hatch)', fontdict={'size':20,'family':'calibri'})
	plt.ylabel('Coefficient of determination ('+"$r^{2}$"+')', fontdict={'size':20,'family':'calibri'})


	plt.style.use('seaborn-white')
	plt.gcf()
	save_name='Top '+ (str)(num_measures)+' ' +control_for+' correlationtrajectories'
	if rescale:
		save_name=save_name+' rescaled'
	plt.savefig(complete_worm_df.save_directory+'/'+save_name+'.pdf')	
	plt.show()
	plt.gcf().clf		
	

def get_all_correlations(measures,first_day,last_day,miRNA,complete_worm_df, save=True,control_for='',rescale=[]):
	
	all_correlations=[]
	all_p_values=[]
	
	for i in range(0, len(measures)):
		correlations,days,p_values,worms,deleted_worms=(get_correlations([measures[i]],first_day,last_day,miRNA,complete_worm_df,control_for,rescale))
		all_correlations.append(correlations)
		all_p_values.append(p_values)
		print("deleted "+(str)(round(100*(len(deleted_worms)/len(worms)),1))+" percent of worms")	
	
	if save:
		dataframe=pandas.DataFrame(data=all_correlations,index=measures,columns=days)
		dataframe.to_csv(complete_worm_df.save_directory+'/correlations.csv')

	return all_correlations,days,all_p_values

def control_for_measure(x_list,y_list,control_means):
	slope_x,intercept_x,rvalue_x,pvalue_x,stderr_x=scipy.stats.linregress(x=control_means,y=y_list)
	residuals=y_list-(control_means*slope_x+intercept_x)
	
	return x_list, residuals


def get_correlations(measure,first_day,last_day,miRNA, complete_worm_df,control_for='',rescale=[]):

	if len(rescale)==0:
		means, days, lifespans, worms, deleted_worms = get_data(complete_worm_df,measure,first_day,last_day)
	else:
		means, days, lifespans, worms, deleted_worms = run_rescaling_multiple(rescale,measure,first_day,last_day)	
	
	if len(control_for)>0 and len(rescale)==0:
		control_means, days, lifespans, worms, deleted_worms=get_data(complete_worm_df, [control_for], first_day, last_day)
	if len(control_for)>0 and len(rescale)>0:
		control_means, days, lifespans, worms, deleted_worms=run_rescaling_multiple(rescale,[control_for],first_day,last_day)
			

	correlations=[]
	p_values=[]
	
	for i in range(0,len(means)):
		if len(control_for)>0:
			x,y=control_for_measure(means[i],lifespans,control_means[i])
			pearson,spearman,yp=run_stats(x,y)
		else:	
			pearson,spearman,yp=run_stats(means[i],lifespans)
		correlations.append(pearson[0]**2)
		p_values.append(pearson[1])

	return correlations,days,p_values,worms,deleted_worms

def plot_all_correlations(first_day,last_day,measures,miRNA,complete_worm_df,show=False, control_for='',rescale=[]):

	for i in range(0,len(measures)):
		plot_correlations(first_day,last_day,[measures[i]],miRNA,complete_worm_df,show,control_for,rescale)


def plot_correlations(first_day,last_day,measure, miRNA, complete_worm_df,show=False,control_for='',rescale=[]):
	correlations,days,p_values,worms,deleted_worms =get_correlations(measure,first_day,last_day,miRNA,complete_worm_df,control_for,rescale)
	print("deleted "+(str)(round(100*(len(deleted_worms)/len(worms)),1))+" percent of worms")	
	plt.style.use('seaborn-white')
	plt.scatter(days,correlations,c=color_terms[measure[0]],marker='o',s=50,edgecolor=color_terms[measure[0]])
	plt.plot(days,correlations, c=color_terms[measure[0]],linestyle='--')
	
	p_x=[]
	p_y=[]
	for i in range(0,len(p_values)):
		if p_values[i]<.05:
			p_y.append(correlations[i]+.03)
			p_x.append(days[i])
	plt.scatter(p_x,p_y,marker=(6,2,0),color=color_terms[measure[0]],s=50)
	if control_for:
		plt.title('Correlation of '+miRNA + '::GFP expression ' + '('+ terms[measure[0]]+')'+ ' with lifespan vs. timepoint, controlled for '+ terms[control_for], y=1.05, fontdict={'size':26,'weight':'bold','family':'calibri'})	
	else:
		plt.title('Correlation of '+miRNA + '::GFP expression ' + '('+ terms[measure[0]]+')'+ ' with lifespan vs. timepoint', y=1.05, fontdict={'size':26,'weight':'bold','family':'calibri'})
	plt.xlabel('Timepoint (day post-hatch)', fontdict={'size':20,'family':'calibri'})
	plt.ylabel('Coefficient of determination ('+"$r^{2}$"+')',fontdict={'size':20,'family':'calibri'})

	ymin, ymax = plt.ylim()
	plt.ylim(0,ymax+.02)
	xmin,xmax=plt.xlim()
	plt.xlim(days[0]-.5,days[-1]+.5)

	plt.gcf()
	if control_for:
		plt.savefig(complete_worm_df.save_directory+'/'+'Correlationtrajectories'+measure[0]+' controlledfor '+ control_for[0]+'.pdf')
	else:	
		plt.savefig(complete_worm_df.save_directory+'/'+'Correlationtrajectories '+measure[0]+'.pdf')	
	if show:
		plt.show(block=False)
		time.sleep(2)
		plt.close()
	plt.gcf().clf		

def make_scatter_plots(all_measures,first_day,last_day, miRNA, complete_worm_df,show=False,control_for='',rescale=[],bin_width_days=2):

	all_means=[]
	
	for i in range(0, len(all_measures)):
		if len(rescale)==0:
			means, days, lifespans, worms, deleted_worms=get_data(complete_worm_df, [all_measures[i]], first_day, last_day)
		if len(rescale)>0:
			means, days, lifespans, worms, deleted_worms=run_rescaling_multiple(rescale, [all_measures[i]], first_day, last_day)	
		if len(control_for)>0:
			control_means, days, lifespans, worms, deleted_worms=get_data(complete_worm_df, [control_for], first_day, last_day)
			residuals_list=[]
			for k in range(0,len(means)):
				x_list,residuals=control_for_measure(means[k],lifespans,control_means[k])
				residuals_list.append(x_list)
		if len(control_for)>0:
			all_means.append(residuals_list)
			lifespans=residuals
		if len(control_for)==0:
			all_means.append(means)
	
	for i in range(0, len(all_means)):
		for j in range(0, len(days)):
			make_scatter_plot(all_means[i][j],lifespans,all_measures[i],days[j],miRNA,complete_worm_df,show,control_for,rescale,bin_width_days=bin_width_days)
	print(deleted_worms)
	print("deleted "+(str)(round(100*(len(deleted_worms)/len(worms)),1))+" percent of worms")		

def make_scatter_plot(means, lifespans, measure, day, miRNA, complete_worm_df,show=False,control_for='',rescale=[],bin_width_days=2):

	life_cohorts, bin_lifes, my_bins, my_colors=sd.life_cohort_bins(complete_worm_df,bin_width_days=bin_width_days)

	if len(rescale)>0:
		my_bins=my_bins/numpy.median(sd.get_lifespans(complete_worm_df)/24)

	pearson,spearman,yp=run_stats(means,lifespans)

	if len(control_for)==0:
		colors=[]
	
		for i in range(0,len(lifespans)):
			for j in range(0,len(my_bins)):
				if lifespans[i]<my_bins[j][1] and lifespans[i]>my_bins[j][0]:
					colors.append(rainbow_colors[j])		
	else:
		colors='grey'				
	plt.scatter(means,lifespans,c=colors, edgecolors=colors,s=35,alpha=.7)
	plt.plot(means,yp, c='gray')
	plt.style.use('seaborn-white')

	title=''
	y_label=''
	save_name=''
	if measure.endswith('gfp') or measure.endswith('rfp'):
		title="Mean "+miRNA+ "::GFP expression at "+ (str)(day) + " dph vs. lifespan"
		x_label="Mean "+miRNA + "::GFP expression ( "+terms[measure]+ " )"
		save_name= "Mean "+miRNA + " expression "+ (str)(day)+ " dph " + measure
		y_label='Lifespan (days)' 	
	if measure.endswith('size'):
		title="Total size at "+ (str)(day) + " dph vs. lifespan"
		x_label="Total size (pixels)"
		save_name= "Total size  "+(str)(day)+ " dph " 
		y_label='Lifespan (days)' 	
	if measure.startswith('intensity'):
		title="Mean autofluorescence at "+ (str)(day) + " dph vs. lifespan"
		x_label="Mean "+terms[measure]
		save_name= "Mean autofluorescence "+(str)(day)+ " dph " + measure
		y_label='Lifespan (days)' 		
	if len(control_for)>0:
		title=title+ ' controlled for '+ terms[control_for]
		save_name=save_name+ 'controlled for '+terms[control_for]
		y_label="Lifespan residual (after regression on "+ terms[control_for] +")"	
	if len(rescale)>0:
		save_name=save_name+' rescaled'
		x_label="Mean "+miRNA + "::GFP expression ( "+terms[measure]+ " )"+", rescaled to median"
		y_label='Lifespan (days post-hatch), rescaled to median' 
	plt.title(title,y=1.05,fontdict={'size':26,'weight':'bold','family':'calibri'})
	plt.ylabel(y_label,fontdict={'size':20,'family':'calibri'})
	plt.xlabel(x_label,fontdict={'size':20,'family':'calibri'})

	if pearson[1]<.00001:
		p="p<.00001"
	else:
		p="p=" + ''+ (str)(round(pearson[1],3))
	if spearman[1]<.00001:
		spearman_p="p<.00001"
	else:
		spearman_p="p=" + '' + (str)(round(spearman[1],3))
				
	ftext="$r^{2}$ = "+(str)(round(pearson[0]**2,3))+" "+p
	gtext=r'$\rho$'+" = "+(str)(round(spearman[0],3))+" "+spearman_p
	more_text="n= " + (str)(len(means))
	
	plt.figtext(.15,.85,ftext,fontsize=20,ha='left')
	plt.figtext(.15,.8,gtext,fontsize=20,ha='left')
	plt.figtext(.8, .2, more_text, fontsize=20, ha='left')
	

	
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

def run_rescaling_multiple(data_frames,measure,first_day,last_day):
	'''		
	Rescale all dataframes individually and then concatenate together
	'''
	all_rescaled_means=[]
	all_lifespans=[]
	all_worms=[]
	all_deleted_worms=[]
	for data_frame in data_frames:
		rescaled_means,days,lifespans,worms,deleted_worms=run_rescaling(data_frame, measure,first_day,last_day)	
		all_rescaled_means.append(rescaled_means)
		all_lifespans.append(lifespans)
		all_worms.append(worms)
		all_deleted_worms.append(deleted_worms)
	all_rescaled_means=numpy.concatenate(all_rescaled_means,axis=1)
	all_lifespans=numpy.concatenate(all_lifespans)
	all_worms=numpy.concatenate(all_worms)
	all_deleted_worms=numpy.concatenate(all_deleted_worms)
	return all_rescaled_means,days,all_lifespans,all_worms,all_deleted_worms

def run_rescaling(complete_worm_df,measure,first_day,last_day):
	'''		
	Rescale a measurement within an individual dataframe to the median of that measurement	
	'''
	means, days, lifespans, worms, deleted_worms=get_data(complete_worm_df,measure,first_day,last_day)
	means=numpy.asarray(means)
	lifespans=numpy.asarray(lifespans)
	rescaled_means=numpy.empty(means.shape)
	rescaled_lifespans=lifespans/numpy.median(lifespans)
	for i in range(0,len(means)):
		median=numpy.median(means[i])
		for j in range(0,len(means[i])):
			rescaled_means[i][j]=means[i][j]/median


	
	return rescaled_means,days,rescaled_lifespans,worms,deleted_worms

def get_data(complete_worm_df,measure,first_day,last_day,rescale=[]):
	#changed select_times on 1/31/18 to be inclusive i.e. 2.0-5.99 dph. Previous version was 2.0 to 6.99 when given 6.0 as last_day argument
	select_times=complete_worm_df.times[complete_worm_df.times.index(first_day):(complete_worm_df.times.index(last_day)+1)]
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
		means[i]=mew.copy()
	
	
	lifespans=numpy.delete(lifespans,nan_positions)
	deleted_worms=[worms[position] for position in nan_positions]
	worms=numpy.delete(worms,nan_positions)

	if len(rescale)>0:
		means=numpy.asarray(means)
		rescaled_means=numpy.empty(means.shape)
		for i in range(0,len(means)):
			median=numpy.median(means[i])
			for j in range(0,len(means[i])):
				rescaled_means[i][j]=means[i][j]/median
		means=rescaled_means.copy()	

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


def measure_slopes(measure, complete_worm_df, first_day, last_day,miRNA,bin_width_days=2):
	worms=complete_worm_df.worms
	lifespans=sd.get_lifespans(complete_worm_df)/24
	select_times=complete_worm_df.times[complete_worm_df.times.index(first_day):(complete_worm_df.times.index(last_day)+1)]
	data=complete_worm_df.mloc(worms,[measure],select_times)
	times=[float(i) for i in select_times]

	life_cohorts, bin_lifes, my_bins, my_colors=sd.life_cohort_bins(complete_worm_df,bin_width_days=bin_width_days)
	slopes=[]

	for i in range(0, len(data)):
		if numpy.any(numpy.isnan(data[i][0])):
			slopes.append(0.00013)	

		else:
			slope,intercept,r,p,std_err=scipy.stats.linregress(times, data[i][0])
			slopes.append(slope)		
	
	slopes=numpy.asarray(slopes)
	lifespans=lifespans[slopes!=0.00013]
	
	print('percent deleted worms is '+(str)(100*(len(worms)-len(lifespans))/len(worms)))

	new_cohorts=[]
	for i in range(0,len(life_cohorts)):
		new_cohorts.append([])
	
	for i in range(0,len(slopes)):
		for j in range(0,len(life_cohorts)):
			if slopes[i]!=0.00013 and i in life_cohorts[j]:
				new_cohorts[j].append(i)
	
	slopes=slopes[slopes!=0.00013]
	new_cohorts=numpy.asarray(new_cohorts)
	
	colors=[]

	for i in range(0,len(worms)):
		for j in range(0,len(new_cohorts)):
			if i in new_cohorts[j]:
				colors.append(rainbow_colors[j])

	plt.style.use('seaborn-white')
	
	plt.title(miRNA +"::GFP slope (days "+(str)(first_day)+" to "+ (str)(last_day)+")" +" vs. lifespan",fontdict={'size':26,'weight':'bold','family':'calibri'})
	plt.xlabel("Slope of "+miRNA+ "::GFP ("+terms[measure]+")",fontdict={'size':20,'family':'calibri'})
	plt.ylabel("Lifespan (days)",fontdict={'size':20,'family':'calibri'})
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

	more_text="n= " + (str)(len(slopes))	
	ftext="$r^{2}$ ="+(str)(round(correlation[0]**2,3))+" "+p
	plt.figtext(.65,.8,ftext,fontsize=20,ha='left')
	plt.figtext(.8, .2, more_text, fontsize=20, ha='left')


	plt.gcf()
	plt.savefig(complete_worm_df.save_directory+'/'+'Slope '+ terms[measure] +(str)(first_day) +' to '+(str)(last_day)+' dph.pdf')	
	plt.show()
	plt.close()
	plt.gcf().clf


def plot_survival_curve(complete_worm_df):
	
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
	plt.plot(days,percent_alive, color='goldenrod')

	plt.title("Survival curve",y=1.05,fontdict={'size':26,'weight':'bold','family':'calibri'})
	plt.xlabel("Lifespan (days)",fontdict={'size':20,'family':'calibri'})
	plt.ylabel("Survival (%)",fontdict={'size':20,'family':'calibri'})
	ftext="mean lifespan="+(str)(round(lifespans.mean(),1)) + " days"
	even_more_text="median lifespan = " + (str)(round(numpy.percentile(lifespans,50),1)) + " days"
	more_text="n= " + (str)(lifespans.size)
	plt.figtext(.5,.8,ftext,fontsize=20,ha='left')
	plt.figtext(.8, .2, more_text, fontsize=20, ha='left')
	plt.figtext(.5, .75, even_more_text, fontsize=20, ha='left')

	plt.gcf()
	plt.savefig(complete_worm_df.save_directory+'/Survival_curve.pdf')
	plt.show()
	plt.gcf().clf

def plot_smoothed_trajectories(data_directories,save_directory, measure, bin_width_days=2, miRNA='',rescale=[]):
	complete_worm_df=make_dataframe(data_directories,save_directory)
	smoothed_data=complete_worm_df.smooth_trajectory(measure,extra_arguments={'smooth_iterations': 3, 'polyorder': 1, 'windowLength': 9, 'smooth_mode': 'savitzky_golay'})
	ages=complete_worm_df.ages
	life_cohorts,bin_lifes,my_bins,my_colors=sd.life_cohort_bins(complete_worm_df,bin_width_days=bin_width_days)
	averages=[]
	
	plt.style.use('seaborn-white')
	plt.xlabel('Age (days)', fontdict={'size':20,'family':'calibri'})

	if measure.endswith('gfp'):
		title='Average '+miRNA+'::GFP profiles of different lifespan cohorts'
		y_label=miRNA+ '::GFP expression ('+ terms[measure]+')'
		save_name= "Average "+terms[measure]+ " profiles"
	if measure.endswith('rfp'):
		title='Average '+miRNA+'::RFP profiles of different lifespan cohorts'
		y_label=miRNA+ '::RFP expression ('+ terms[measure]+')'
		save_name= "Average "+terms[measure]+ " RFP profiles"	
	if measure.endswith('size'):
		title="Average size of different lifespan cohorts"
		y_label="Size"
		save_name= "Average size profiles" 	
	if measure.startswith('intensity'):
		title="Average autofluorescence profiles of different lifespan cohorts"
		y_label="Autofluorescence ("+measure+")" 	
		save_name= "Average "+terms[measure]+" profiles"	

	plt.title(title,y=1.05,fontdict={'size':26, 'weight':'bold','family':'calibri'})
	plt.ylabel(y_label, fontdict={'size':20,'family':'calibri'})
	more_text='n =  '
	
	for i in range(0, len(life_cohorts)):
		more_text=more_text+'       '+ (str)(len(life_cohorts[i]))

	plt.figtext(.5, .2, more_text, fontsize=20, ha='left')
	
	for i in range(0, len(life_cohorts)):
		cohort_data = smoothed_data[life_cohorts[i]]
		cohort_data = cohort_data[~numpy.isnan(cohort_data).all(axis = 1)]
		data=numpy.mean(cohort_data, axis = 0)
		mew=numpy.copy(data)
		if numpy.searchsorted(ages,bin_lifes[i])<len(ages):
			mew[(numpy.searchsorted(ages,bin_lifes[i])):]='NaN'
		#mew[ages.index(round(bin_lifes[i]),1):]='NaN'
		averages.append(mew)
		
	for i in range(0, len(life_cohorts)):
		plt.plot(ages, averages[i], color=rainbow_colors[i])
	ymin,ymax=plt.ylim()
	plt.ylim(0,ymax+.1*ymax)	

	
	plt.gcf()
	plt.savefig(complete_worm_df.save_directory+'/'+save_name+' smoothed.pdf')
	plt.show(block=False)
	time.sleep(2)
	plt.close()
	plt.gcf().clf
	return ages, averages
def plot_all_trajectories(data_directories,save_directory, measures, bin_width_days=2,miRNA='',smooth=False):
	for i in range(0,len(measures)):
		if smooth:
			plot_smoothed_trajectories(data_directories,save_directory, measures[i],bin_width_days,miRNA)
		else:		
			plot_trajectories(data_directories,save_directory, measures[i],bin_width_days, miRNA)	

def plot_trajectories(data_directories,save_directory, measure,bin_width_days=2, miRNA=''):
	complete_worm_df=make_dataframe(data_directories,save_directory)
	ages=complete_worm_df.ages
	life_cohorts,bin_lifes,my_bins,my_colors=sd.life_cohort_bins(complete_worm_df,bin_width_days=bin_width_days)
	averages=[]
	
	if measure.endswith('gfp'):
		title='Average '+miRNA+'::GFP profiles of different lifespan cohorts'
		y_label=miRNA+ '::GFP expression ('+ terms[measure]+')'
		save_name= "Average "+terms[measure]+ " profiles"
	if measure.endswith('size'):
		title="Average size of different lifespan cohorts"
		y_label="Size"
		save_name= "Average size profiles" 	
	if measure.startswith('intensity'):
		title="Average autofluorescence profiles of different lifespan cohorts"
		y_label="Autofluorescence ("+measure+")" 	
		save_name= "Average "+terms[measure]+" profiles"	
	
	plt.style.use('seaborn-white')
	plt.xlabel('Age (days)', fontdict={'size':20,'family':'calibri'})
	
	plt.title(title,y=1.05,fontdict={'size':26, 'weight':'bold','family':'calibri'})
	plt.ylabel(y_label, fontdict={'size':20,'family':'calibri'})
	more_text='n =  '

	for i in range(0, len(life_cohorts)):
		more_text=more_text+'       '+ (str)(len(life_cohorts[i]))

	plt.figtext(.5, .2, more_text, fontsize=20, ha='left')
	
	for i in range(0, len(life_cohorts)):
		averages.append(sd.cohort_trace(complete_worm_df,life_cohorts[i],measure,False))
	
	for i in range(0, len(life_cohorts)):
		
		plt.plot(ages, averages[i], color=rainbow_colors[i])
		
	save_name='Average'+measure+'profiles'	

	plt.gcf()
	plt.savefig(complete_worm_df.save_directory+'/'+save_name+'.pdf')
	plt.show(block=False)
	time.sleep(1)
	plt.close()
	plt.gcf().clf		
	return ages, averages





