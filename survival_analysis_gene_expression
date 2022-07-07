for i in  list(genes_apopt.columns):  #list gene  

    # convert value of gene expression in 0,1,2 using quantile method
    data_surv_zscore.loc[data_surv_zscore[i] > float(data_surv_zscore[i].quantile([0.75])),i] = 2 # > quartile 0.75 

    data_surv_zscore.loc[((data_surv_zscore[i] > float(data_surv_zscore[i].quantile([0.25]))) & #0.25 to 0.75 quartile
                     (data_surv_zscore[i] < float(data_surv_zscore[i].quantile([0.75])))),i]=1

    data_surv_zscore.loc[data_surv_zscore[i] < float(data_surv_zscore[i].quantile([0.25])),i] = 0 # < 0.25 
    # Fit Kelpen Meier
    kmf_high = KaplanMeierFitter()
    kmf_med = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()
    
    #create variable samples with high, med and low 
    High_gene=data_surv_zscore.query(f"{i} == 2")
    Med_gene=data_surv_zscore.query(f"{i} == 1")
    Low_gene=data_surv_zscore.query(f"{i} == 0")
    
    #Fitting
    kmf_high.fit(durations=High_gene["Event.Free.Survival.Time.in.Days"],event_observed=High_HES4["Dead"], label=f"{i}_high")
    kmf_med.fit(durations=Med_gene["Event.Free.Survival.Time.in.Days"],event_observed=Med_HES4["Dead"], label=f"{i}_med")
    kmf_low.fit(durations=Low_gene["Event.Free.Survival.Time.in.Days"],event_observed=Low_HES4["Dead"], label=f"{i}_low")
       
    #High vs MED #stats significative
    print(f"High vs Med of {i} gene")
    T=High_gene["Event.Free.Survival.Time.in.Days"]
    E=High_gene["Dead"]
    T1=Med_gene["Event.Free.Survival.Time.in.Days"]
    E1=Med_gene["Dead"]
    results= logrank_test(T,T1,event_observed_A=E,event_observed_B=E1)
    results.print_summary()
    
    #High vs Low #stats significative
    print(f"\n\n High vs Low of {i} gene")
    T=High_gene["Event.Free.Survival.Time.in.Days"]
    E=High_gene["Dead"]
    T1=Low_gene["Event.Free.Survival.Time.in.Days"]
    E1=Low_gene["Dead"]
    results= logrank_test(T,T1,event_observed_A=E,event_observed_B=E1)
    results.print_summary()
    
    #Med vs Low
    print(f"\n\n Med vs Low of {i} gene")
    T=Med_gene["Event.Free.Survival.Time.in.Days"]
    E=Med_gene["Dead"]
    T1=Low_gene["Event.Free.Survival.Time.in.Days"]
    E1=Low_gene["Dead"]
    results= logrank_test(T,T1,event_observed_A=E,event_observed_B=E1)
    results.print_summary()
    
    #create and show plot
    kmf_high.plot()
    kmf_med.plot()
    kmf_low.plot()
    plt.gcf().set_dpi(150)
    plt.xlabel("Days passed")
    plt.ylabel("Survival Probability")
    plt.title("KMF")
    plt.show()
