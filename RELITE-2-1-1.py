                #++++++++++++++++++++++++++++++++
                #       Welcome to RELITE       #
                #++++++++++++++++++++++++++++++++
# To obtain netcdf files of output graphs replace respective 
# @ct.output.figure() or livefigure() from line 44-63 with 
# @ct.output.download(). 
# Also post processing scripts and aditional data for figure 10 and 11 in
# publication will be made available by request via email.
# My email address : soineade@gmail.com

import cdstoolbox as ct

app_layout = {
    'input_ncols' : 2,
    'output_ncols' : 2,
    'output_align': 'bottom',
}

@ct.application(title='Renewable Electricity Synergy (RELITE): Solar & Wind',layout=app_layout)

@ct.input.text('area',label='Domain',default='90/-180/-90/180', description= 'Enter coordinates:[North/West/South/East]')

@ct.input.text('grid',label='Grid Size',default='3/3',
 description= 'Enter grid size')

@ct.input.text('device',label='Wind Turbine and Solar Panel Characteristics',default='117/3/12/22.5',
 description= 'Enter Wind turbine properties (Hub height/Cut-in speed/rated wind speed/cut-out speed). Note: The solar PV CF is modelled for only monocrystalline silicon cell')

@ct.input.text( 'lat_lon',label='Select S-W hybrid power station ',default='-27/117',
description='Enter latitude/Longitude in decimal degrees')

@ct.input.slider('year', min = 1979, max = 2021, step = 1,default=(2011,2011),label='Time Range', description='Note: Time range maximum (upper boundary) should be changed  when the maximum time becomes obsolate. . This can be done in the input.slider section of the code')

@ct.input.text('utc',label='UTC to Local time ',default='8', description= 'Enter to shift Utc by defined value')

@ct.input.text('ratio',label='Capacity ratio (n:m) for Stability coefficient ',default='1:1',
 description= 'Enter desired capacity ratio (power mix) for solar(n) and wind power(m) to enable calculation Cstab')

@ct.input.text('thres',label='Synergy Threshold',default='0.4/0.2/0.7/0.5',
 description= 'Enter threshold for good dinurnal (cstab) and seasonal (Pearson) synergy in the format cstab/pearson')

@ct.input.checkbox('CF',['wind capacity factor','solar capacity factor'],label=' Average capacity factor map', link={'wind capacity factor': 'output-0', 'solar capacity factor':'output-1',})

@ct.input.checkbox('Metric',['Stability coefficient (hourly)','Pearson correlation coefficient (monthly)'],label='Metric (Spatial synergy)', link={'Stability coefficient (hourly)': 'output-2','Pearson correlation coefficient (monthly)': 'output-3'},description= 'Note: For seasoanal synergy Pearson correlation  is used.')


@ct.input.checkbox('good_sy',['Synergy Zones',],label='Qualitavtive Map', link={'Synergy Zones': 'output-4',},description= 'Bad synergic zones are masked out')

@ct.input.checkbox('cstab_mon',['Monthly Average Stability Coefficicent','Monthly power profile'],label=' Temporal Synergy', link={'Monthly Average Stability Coefficicent': 'output-5','Monthly power profile':'output-6'},description= 'Produces a line graph of monthly averaged Cstab and power profile for defined latitude and longitude coordinate')

@ct.output.figure()
@ct.output.figure()
@ct.output.figure()
@ct.output.figure()
@ct.output.figure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()
@ct.output.livefigure()


def application(area,utc,grid,device,CF,Metric,year,lat_lon,cstab_mon,thres,good_sy,ratio):
       
    #====================================
    #Extrating geogprahic extent 
    #====================================
    #print(area)
    print(lat_lon)
    extent=[]
    for i in area.split('/'):
        extent.append(int(i))
    print(extent)

    #=========================================
    #Extrating wind turbine properties as float
    #========================================= 
    windprop=[]
    for i in device.split('/'):
        windprop.append(float(i))
    #=========================================
    #year length
    #=========================================
    Length=int(year[1])-int(year[0])
    #start
    start_year=int(year[0])
    #end (always+1)
    end_year=int(year[0])+Length+1
    #loop for year list
    List_year=[]
    for i in range(start_year,end_year):
             List_year.append(i)          
    year=List_year
    print(year)
    #==========================================
    Cstab_all=[]
    CF_solar=[]
    CF_wind=[]
    for i in range(len(year)):
            print(year[i])
            U100,V100,U10,V10,T = ct.catalogue.retrieve(
                    'reanalysis-era5-single-levels',
                    {
                        'product_type': 'reanalysis',
                        'variable': [
                            '100m_u_component_of_wind', '100m_v_component_of_wind', '10m_u_component_of_wind',
                            '10m_v_component_of_neutral_wind', '2m_temperature'],
                        'year': year[i], 
                        'month': [
                            '01', '02', '03',
                            '04', '05', '06',
                            '07', '08', '09',
                            '10', '11', '12',
                        ],
                        'day': [
                            '01', '02', '03',
                            '04', '05', '06',
                            '07', '08', '09',
                            '10', '11', '12',
                            '13', '14', '15',
                            '16', '17', '18',
                            '19', '20', '21',
                            '22', '23', '24',
                            '25', '26', '27',
                            '28', '29', '30',
                            '31',
                        ],
                        'time': [
                            '00:00', '01:00', '02:00',
                            '03:00', '04:00', '05:00',
                            '06:00', '07:00', '08:00',
                            '09:00', '10:00', '11:00',
                            '12:00', '13:00', '14:00',
                            '15:00', '16:00', '17:00',
                            '18:00', '19:00', '20:00',
                            '21:00', '22:00', '23:00',
                        ],
                        'area': extent,
                        'grid': grid,
                    }
                )            
            #+++++++++++++++++++++++++++++++++++++++++++++++++++
            #New Additional format for SSRD variable (on DEC 07)
            #+++++++++++++++++++++++++++++++++++++++++++++++++++
            S = ct.catalogue.retrieve(
                    'reanalysis-era5-single-levels',
                    {
                        'product_type': 'reanalysis',
                        'variable':'surface_solar_radiation_downwards',
                        'year': year[i], 
                        'month': [
                            '01', '02', '03',
                            '04', '05', '06',
                            '07', '08', '09',
                            '10', '11', '12',
                        ],
                        'day': [
                            '01', '02', '03',
                            '04', '05', '06',
                            '07', '08', '09',
                            '10', '11', '12',
                            '13', '14', '15',
                            '16', '17', '18',
                            '19', '20', '21',
                            '22', '23', '24',
                            '25', '26', '27',
                            '28', '29', '30',
                            '31',
                        ],
                        'time': [
                            '00:00', '01:00', '02:00',
                            '03:00', '04:00', '05:00',
                            '06:00', '07:00', '08:00',
                            '09:00', '10:00', '11:00',
                            '12:00', '13:00', '14:00',
                            '15:00', '16:00', '17:00',
                            '18:00', '19:00', '20:00',
                            '21:00', '22:00', '23:00',
                        ],
                        '_rate_type' : 'rate',
                        'area': extent,
                        'grid': grid,
                    }
                )
            
            
            T = ct.cdm.update_attributes(T, remove=['standard_name','units'])
            T=ct.operator.sub(T,273.15)#convert to degrees celcius
            S = ct.cdm.update_attributes(S, remove=['standard_name','units'])
            #==================================================
            #Calculating the net velocity at height 10 and 100m
            #==================================================
            Vnet10=ct.math.sqrt(ct.operator.add(ct.math.square(U10),ct.math.square(V10)))
            Vnet100=ct.math.sqrt(ct.operator.add(ct.math.square(U100),ct.math.square(V100)))

            #=====================================================
            # Estimating Rougness height(Hr)function for Each hour 
            #=====================================================
            h1=100
            h0=10
            b=ct.operator.sub(Vnet10,Vnet10)
            h1=ct.cube.where(ct.math.logical_not(b<=0),b,h1)
            h0=ct.cube.where(ct.math.logical_not(b<=0),b,h0)
            v0=Vnet10
            v1=Vnet100
            Hr=ct.math.exp(ct.operator.div(ct.operator.sub(ct.operator.mul(v0,ct.math.log(h1)),ct.operator.mul(v1,ct.math.log(h0))),(ct.operator.sub(v0,v1))))
            HrNew=ct.math.nan_to_num(Hr)
            #===================================================
            ##Estimating wind velocity  at hub height 117m (V117)
            #===================================================
            h1=windprop[0]#hub height 
            h0=100
            h1=ct.cube.where(ct.math.logical_not(b<=0),b,h1)
            h0=ct.cube.where(ct.math.logical_not(b<=0),b,h0)

            V117_N=ct.operator.div(ct.math.log(ct.operator.div(h1,HrNew))
            ,ct.math.log(ct.operator.div(h0,HrNew)))
            
            V117=ct.operator.mul(Vnet100,V117_N)
            V117=ct.math.nan_to_num(V117)
            V117=ct.cdm.update_attributes( V117, remove=['standard_name','units'])

            #============================================================
            #Calulating wind Capacity Factor (WCF) Note:windspeed in ms-1
            #============================================================
            Vin=windprop[1]#cut in wind speed 
            Vr=windprop[2]#rated wind speed
            Vc=windprop[3]#cut off wind speed
            
            #Data conditioning 
            WCF=ct.cube.where(ct.math.logical_not(V117<=Vin),V117,0)  
            WCF=ct.cube.where(ct.math.logical_not(ct.math.logical_and(V117>Vin,V117<=Vr)),WCF,ct.operator.div(ct.operator.sub(ct.operator.pow(V117,3),Vin**3),
(Vr**3-Vin**3))) 
            WCF=ct.cube.where(ct.math.logical_not(ct.math.logical_and(V117>Vr,V117<=Vc)),WCF,1) 
            WCF=ct.cube.where(ct.math.logical_not(V117>Vc),WCF,0)
            WCF=ct.cube.where(ct.math.logical_not(WCF>1),WCF,1)#
            WCF=ct.cube.where(ct.math.logical_not(WCF<0),WCF,0)
            #================================================
            #Wind Capacity Factor data for synergy calculation 
            #================================================    
            CFw=WCF # will be used for synergy calculation
            #==========================================+=====
            
            CF_wind.append(CFw)
                                       
            #SOLAR CAPACITY FACTOR(SCF)*** 
            #==================
            #General parameters
            #==================
            Beta=0.0045
            gamma=0.1           
            Tref=25#◦C
            Gref=1000
            c1=-3.75
            c2=1.14
            c3=0.0175

            
            #====================================================
            #calculating hourly SCF  
            #====================================================
            Tcell=c1+c2*T+c3*S
            Sr=ct.math.log10(S)
            Sr=ct.math.nan_to_num(Sr)
            Sr=ct.cube.where(ct.math.logical_not(Sr<0),Sr,0)
            Sr= ct.cdm.update_attributes(Sr, remove=['standard_name','units'])           
            g_sr=gamma*Sr
            Ncell_B=Beta*(Tcell-Tref)
            Ncell=1-Ncell_B+g_sr
            Ncell=ct.cube.where(ct.math.logical_not(Ncell>1),Ncell,1)
            Ncell= ct.cdm.update_attributes( Ncell, attrs={'units': ''})            
            SCF=(ct.operator.mul( Ncell,S))/(Gref)
            SCF=ct.cube.where(ct.math.logical_not(SCF>1),SCF,1)#
            SCF=ct.cube.where(ct.math.logical_not(SCF<0),SCF,0)
            
            
            #==================================================
            #Solar Capacity Factor data for synergy calculation 
            #==================================================
            CFs=SCF #will be used for synergy calculation
            #CFs=ct.cdm.update_attributes(CFs, attrs={'long_name': ''})
            #CFs=ct.cdm.update_attributes(CFs, remove=['standard_name','units']) 
            #==================================================
            
            #SCF=ct.cube.average(SCF, dim='time')
            CF_solar.append(CFs)
            
            #==================================================
            #Calculating cstab daily
            #==================================================
            n_m=[]
            for i in ratio.split(':'):
                n_m.append(float(i))
            print(n_m)
            n=n_m[0]
            m=n_m[1]
            CFadd_t=ct.operator.add(ct.operator.mul(n,CFs),ct.operator.mul(m,CFw))
            CFmix_t=ct.operator.div(CFadd_t,(n+m))
            
            CFmix_Avg=ct.cube.resample(CFmix_t, freq='day', dim='time',how='mean')
            CFs_daily_Avg=ct.cube.resample(CFs, freq='day', dim='time',how='mean')

            #CFmix_anomaly
            sel=CFmix_t
            dailymean=ct.climate.climatology_mean(sel,frequency='dayofyear')
            sub= ct.climate.anomaly(sel, climatology= dailymean,frequency='dayofyear')
            CSub=sub            
            #========
            #SSQ calculation
            #==========
            #Numerator
            Csubsq=ct.math.square(CSub)
            Num_sum=ct.cube.resample(Csubsq, freq='day',dim='time',how='sum')

            Num_sum_sqrt=ct.math.sqrt(Num_sum)
           
           #=============
           #denominator
           #===========
            #anomaly_CFs
            sel=CFs
            dailymean= ct.climate.climatology_mean(sel, frequency='dayofyear')
            sub_CFs= ct.climate.anomaly(sel, climatology= dailymean, frequency='dayofyear')
            CSubden=sub_CFs
            #================
            #SSQ calculation
            #===============
            Csubsqden=ct.math.square(CSubden)
            Den_sum=ct.cube.resample (Csubsqden, freq='day', dim='time', how='sum')
            Den_sum_sqrt=ct.math.sqrt(Den_sum)###
            #=====================
            #Dividing Num and den
            #=====================
            CFdiv_ssq=ct.operator.div(Num_sum_sqrt,Den_sum_sqrt)
            CFdiv=ct.operator.div(CFs_daily_Avg,CFmix_Avg)
            Cstab_div=ct.operator.mul(CFdiv_ssq,CFdiv)
            Cstab=ct.operator.sub(1,Cstab_div)         
            Cstab_all.append(Cstab)
   
     
    #==============================================
    #concatenating all CFs to full years
    #==============================================
    for i in range(len(CF_solar)-1):
                   CF_solar[i+1] = ct.cube.concat([CF_solar[i],CF_solar[i+1]],dim='time')
    SCF_full=CF_solar[-1]
    #taking avreage over time dimension 
    SCF=ct.cube.average(SCF_full,dim='time')
    SCF=ct.cdm.update_attributes(SCF, remove=['standard_name','units'])
    SCF=ct.cdm.update_attributes(SCF,attrs={'units':''})
    SCF=ct.cdm.update_attributes(SCF, attrs={'long_name':'CF'})
    
    #==============================================
    #concatenating all CFw to full years
    #==============================================
    for i in range(len(CF_wind)-1):
                   CF_wind[i+1] = ct.cube.concat([CF_wind[i],CF_wind[i+1]],dim='time')
    WCF_full=CF_wind[-1]
 
    #taking avreage over time dimension 
    WCF=ct.cube.average(WCF_full,dim='time')
    WCF=ct.cdm.update_attributes(WCF, remove=['standard_name','units'])
    WCF=ct.cdm.update_attributes(WCF,attrs={'units':''})
    WCF=ct.cdm.update_attributes(WCF, attrs={'long_name':'CF'})
    
    p_area=ct.geo.cell_area(SCF)
    
    #=========================================================       
    #concatenating all cstab to full years
    #=========================================================       
    for i in range(len(Cstab_all)-1):
                   Cstab_all[i+1] = ct.cube.concat([Cstab_all[i],Cstab_all[i+1]],dim='time')
    Cstab_full=Cstab_all[-1]
    #taking avreage over time dimension 
    Cstab_Avg=ct.cube.average(Cstab_full,dim='time')
    Cstab_Avg=ct.cdm.update_attributes(Cstab_Avg, remove=['standard_name','units'])
    Cstab_Avg=ct.cdm.update_attributes(Cstab_Avg,attrs={'units': ''})
    Cstab_Avg=ct.cdm.update_attributes(Cstab_Avg, attrs={'long_name': 'Cstab'})
    Cstab_Avg=ct.cube.select(Cstab_Avg, extent=[-180, 180, -60, 60])
        
    print('max_Scf_'+str(ct.math.nanmax(SCF)))
    print('max_wcf_'+str(ct.math.nanmax(WCF)))

    #==================================================================
    #                         Maps 
    # Note !!! oceans are only shaded as white  but also contains 
    # values (required for offshore analysis of solar and wind power)
    #==================================================================

    MAP_CONFIG = {
        'title': {
        'text_lines':['Wind Capacity Factor (-)'],
        'text_justification': 'centre',
        'text_font_size': 0.85,
        'text_colour': 'charcoal'},
        'projection':{
    'subpage_map_projection': 'cylindrical',
    'subpage_lower_left_latitude': extent[2],
    'subpage_lower_left_longitude':extent[1],
    'subpage_upper_right_latitude':extent[0],
    'subpage_upper_right_longitude': extent[3]
},
    'contour': {
        'legend':'on',
        'contour': 'off',
        'contour_shade': 'on',
        'contour_shade_method': 'area_fill',
        'contour_shade_colour_method': 'palette',
        'contour_shade_palette_name': 'eccharts_yellow_red_9',
        'contour_level_selection_type': 'count',
        'contour_level_count': 10,
        #fix magics lower and upper baoundary interpolation error
        'contour_interpolation_floor':0,
        'contour_interpolation_ceiling':ct.math.nanmax(WCF),
        'contour_min_level':0,
        'contour_max_level':'off',
        'contour_label': 'off',

        
    },
    'background': {
        'map_grid' : "off",
        'map_coastline_sea_shade': 'on',
        'map_coastline_sea_shade_colour': 'white'
    },
        'foreground': {
            
        'map_grid' : "off",
        'map_coastline_sea_shade': "on",
        'map_coastline_sea_shade_colour': "white",
        'map_coastline_land_shade': "off",
        'map_coastline_land_shade_colour': "cream",
        'map_boundaries': "on",
        'map_boundaries_colour': "black",
        'map_label_height': 0.4
    },
    
    'legend': {
        'legend_text_colour': "black",
        'legend_text_font_size': '50%',
        'legend_display_type':"continuous",
    },
    'width': 900
}

    #windCF
    fig_WCF = ct.map.plot(WCF,**MAP_CONFIG)
     
    #update solar CF title and plot 
    MAP_CONFIG['title']['text_lines']=["Solar Capacity Factor (-)"]
    MAP_CONFIG['contour'][ 'contour_interpolation_ceiling']=ct.math.nanmax(SCF)
    fig_SCF =  ct.map.plot(SCF,**MAP_CONFIG)    
    
    
    #update cstab title and attributes 
    MAP_CONFIG['title']['text_lines']=["Diurnal Synergy (-)"]
    MAP_CONFIG['contour'][ 'contour_interpolation_ceiling']=ct.math.nanmax(Cstab_Avg)
    MAP_CONFIG['contour']['contour_max_level']=ct.math.nanmax(Cstab_Avg)
    
        
    fig_Cstab =ct.map.plot(Cstab_Avg,**MAP_CONFIG)
    
    
    #=================================================
    #Obtaining  lon and lat for power profile and cstab
    #==================================================
    Lat_lon=[]
    for i in lat_lon.split('/'):
        Lat_lon.append(float(i))
    lat=Lat_lon[0]
    lon=Lat_lon[1]
     
    #==================================================================
    #             Annual average solar & wind CF
    # Please note that annual average values are printed in the console
    #==================================================================
    print('Annual average solar & wind CF')
    
    #Mean SCF
    an_avg_solar=ct.geo.extract_point( SCF, lon=lon, lat=lat)
    an_avg_solar=ct.cdm.get_value(an_avg_solar)
    print('Annual average solar CF='+str(an_avg_solar))
    #Mean WCF
    an_avg_wind=ct.geo.extract_point( WCF, lon=lon, lat=lat)
    an_avg_wind=ct.cdm.get_value(an_avg_wind)
    print('Annual average wind CF='+str(an_avg_wind))
    
    #==================================================================
    #     Annual average capacity factors for standard deviation 
    # Please note that annual average capacity factors are printed to   
    # the console which are then used to calculate the standard deviation
    # of capicity factors offline 
    #==================================================================
    # standard deviation SCF
    SCF_an_full=[]
    for i in range(len(year)):
        SCF_an=ct.cube.select(SCF_full, time=str(year[i]))
        SCF_an_av=ct.cube.average(SCF_an, dim='time')
        SCF_forstd=ct.geo.extract_point( SCF_an_av , lon=lon, lat=lat)
        SCF_an_full.append(SCF_forstd)
    print('Annual average data for STD of solar CF is')
    print(SCF_an_full)

    #standard deviation WCF
    WCF_an_full=[]
    for i in range(len(year)):
        WCF_an=ct.cube.select(WCF_full, time=str(year[i]))
        WCF_an_av=ct.cube.average(WCF_an, dim='time')
        WCF_forstd=ct.geo.extract_point(WCF_an_av, lon=lon, lat=lat)
        WCF_an_full.append(WCF_forstd)
    print('Annual average data for STD  of wind CF is')
    print(WCF_an_full)
  
    #==============================================
    #Plots monthly average cstab 
    #===============================================
    Cstab_month=ct.cube.groupby_reduce(Cstab_full, group='time.month', how='mean', dim='time')
    Cstab_mon=ct.geo.extract_point( Cstab_month, lon=lon, lat=lat)
    Cstab_error=ct.geo.extract_point(Cstab_full, lon=lon, lat=lat) 
    Cstab_month_error=ct.climate.monthly_std(Cstab_error)
    
    #==============================================
    #Plots monthly average cstab 
    #==============================================
    xaxis={'title':'Month','tickvals':[1,2,3,4,5,6,7,8,9,10,11,12],'ticktext':['J','F','M','A','M','J','J','A','S','O','N','D'],'ticks':'outside','showline':True,'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18}
    
    Cstab_ma= ct.chart.line(Cstab_mon,error_y=  Cstab_month_error,
    scatter_kwargs={'name':'Cstab'
    },layout_kwargs = {
         'title': ' Monthly average cstab',
         'xaxis':xaxis,
         'yaxis': {'title':'Cstab(-)','range': [0,1],             'ticks':'outside','showline':True,'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18}})
    
    #====================================
    # Solar:  hourly average CF by month
    #=====================================
    selcfs=ct.geo.extract_point(SCF_full, lon=lon, lat=lat)
    selcfs=ct.cube.shift_coordinates(selcfs, {'time': utc + 'hours'})
    month_cfs=[]
    print('hourly average_solar 50%ile')
    for i in range(12):
        hr_s_av=ct.climate.season_select(selcfs,rule='month',start=i+1, stop=i+1)       
        hourly_average_solar = ct.cube.groupby_reduce(
        data=hr_s_av,
        group='time.hour',
        how='quantile',q=0.5,
        dim='time')
        month_cfs.append(hourly_average_solar)
        
    month_cfs_25=[]
    print('hourly average_solar_25%tile')
    for i in range(12):
        hr_s_av=ct.climate.season_select(selcfs,rule='month',start=i+1, stop=i+1)       
        hourly_average_solar = ct.cube.groupby_reduce(
        data=hr_s_av,
        group='time.hour',
        how='quantile',q=0.25,
        dim='time')
        month_cfs_25.append(hourly_average_solar)
        
    month_cfs_75=[]
    print('hourly average_solar_75%tile')
    for i in range(12):
        hr_s_av=ct.climate.season_select(selcfs,rule='month',start=i+1, stop=i+1)       
        hourly_average_solar = ct.cube.groupby_reduce(
        data=hr_s_av,
        group='time.hour',
        how='quantile',q=0.75,
        dim='time')
        month_cfs_75.append(hourly_average_solar)
    
    #===============================================
    # Wind: hourly average capacity factors by month
    #================================================
    selwfs=ct.geo.extract_point(WCF_full, lon=lon, lat=lat)
    selwfs=ct.cube.shift_coordinates(selwfs, {'time': utc + 'hours'})
    month_wfs=[]
    print('hourly average_wind')
    for i in range(12):
        hr_w_av= ct.climate.season_select(selwfs,rule='month',start=i+1, stop=i+1)
        hourly_wind = ct.cube.groupby_reduce(
        data= hr_w_av,
        group='time.hour',
        how='quantile',q=0.5,
        dim='time')   
        month_wfs.append(hourly_wind)
      
    month_wfs_25=[]
    print('hourly average_wind_25%tile')
    for i in range(12):
        hr_w_av= ct.climate.season_select(selwfs,rule='month',start=i+1, stop=i+1)
        hourly_wind = ct.cube.groupby_reduce(
        data= hr_w_av,
        group='time.hour',
        how='quantile',q=0.25,
        dim='time')   
        month_wfs_25.append(hourly_wind)
    
    month_wfs_75=[]
    print('hourly average_wind_75%tile')
    for i in range(12):
        hr_w_av= ct.climate.season_select(selwfs,rule='month',start=i+1, stop=i+1)
        hourly_wind = ct.cube.groupby_reduce(
        data= hr_w_av,
        group='time.hour',
        how='quantile',q=0.75,
        dim='time')   
        month_wfs_75.append(hourly_wind)
    #===================================================================
    # Plots  for houly average by month (solar and wind capacity factors)
    #====================================================================
    wind={'color':'blue'}
    wind_q={'color':'rgba(0, 0, 255, 0.2)'}
    solar={'color':'red'}
    solar_q={'color':'rgba(255, 0, 0, 0.2)'}
    xaxis={'title':'Hour(UTC)','tickvals':[0,12,23],'ticktext':['0', '12', '23',], 'ticks':'outside','showline':True,
       'linewidth':1,'linecolor':'black','mirror':True, 'tickfont_size':18}
            
    #==================================================
    #Jananry 
    #==================================================
    jan= ct.chart.line( month_cfs[0],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'Jan','yaxis': {'title':'CF(-)','range': [0,1],'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    jan= ct.chart.line( month_cfs_75[0],line=solar_q,fig= jan,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    jan= ct.chart.line( month_cfs_25[0],line=solar_q,fig=jan,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    jan= ct.chart.line( month_wfs[0],line=wind,fig= jan,scatter_kwargs={'mode':'lines','name':'wind' })   
    jan= ct.chart.line( month_wfs_75[0],line=wind_q,fig= jan,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    jan= ct.chart.line( month_wfs_25[0],line=wind_q,fig= jan,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
    
    #==================================================
    #February 
    #==================================================
    feb= ct.chart.line( month_cfs[1],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'Feb','yaxis': {'title':'CF(-)','range': [0,1], 'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    feb= ct.chart.line( month_cfs_75[1],line=solar_q,fig= feb,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    feb= ct.chart.line( month_cfs_25[1],line=solar_q,fig=feb,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    feb= ct.chart.line( month_wfs[1],line=wind,fig= feb,scatter_kwargs={'mode':'lines','name':'wind' })   
    feb= ct.chart.line( month_wfs_75[1],line=wind_q,fig= feb,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    feb= ct.chart.line( month_wfs_25[1],line=wind_q,fig= feb,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
    
    #==================================================
    #March
    #==================================================
    mar= ct.chart.line( month_cfs[2],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'Mar','yaxis': {'title':'CF(-)','range': [0,1], 'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    mar= ct.chart.line( month_cfs_75[2],line=solar_q,fig= mar,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    mar= ct.chart.line( month_cfs_25[2],line=solar_q,fig=mar,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    mar= ct.chart.line( month_wfs[2],line=wind,fig= mar,scatter_kwargs={'mode':'lines','name':'wind' })   
    mar= ct.chart.line( month_wfs_75[2],line=wind_q,fig= mar,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    mar= ct.chart.line( month_wfs_25[2],line=wind_q,fig= mar,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
     
    #==================================================
    #April
    #==================================================
    apr= ct.chart.line( month_cfs[3],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'Apr','yaxis': {'title':'CF(-)','range': [0,1],'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    apr= ct.chart.line( month_cfs_75[3],line=solar_q,fig= apr,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    apr= ct.chart.line( month_cfs_25[3],line=solar_q,fig=apr,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    apr= ct.chart.line( month_wfs[3],line=wind,fig= apr,scatter_kwargs={'mode':'lines','name':'wind' })   
    apr= ct.chart.line( month_wfs_75[3],line=wind_q,fig= apr,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    apr= ct.chart.line( month_wfs_25[3],line=wind_q,fig= apr,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
    
    #==================================================
    #May
    #==================================================
    may= ct.chart.line( month_cfs[4],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'May','yaxis': {'title':'CF(-)','range': [0,1],'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    may= ct.chart.line( month_cfs_75[4],line=solar_q,fig= may,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    may= ct.chart.line( month_cfs_25[4],line=solar_q,fig=may,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    may= ct.chart.line( month_wfs[4],line=wind,fig= may,scatter_kwargs={'mode':'lines','name':'wind' })   
    may= ct.chart.line( month_wfs_75[4],line=wind_q,fig= may,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    may= ct.chart.line( month_wfs_25[4],line=wind_q,fig= may,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
    
    #==================================================
    #June
    #==================================================
    jun= ct.chart.line( month_cfs[5],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'Jun','yaxis': {'title':'CF(-)','range': [0,1],'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    jun= ct.chart.line( month_cfs_75[5],line=solar_q,fig= jun,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    jun= ct.chart.line( month_cfs_25[5],line=solar_q,fig=jun,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    jun= ct.chart.line( month_wfs[5],line=wind,fig= jun,scatter_kwargs={'mode':'lines','name':'wind' })   
    jun= ct.chart.line( month_wfs_75[5],line=wind_q,fig= jun,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    jun= ct.chart.line( month_wfs_25[5],line=wind_q,fig= jun,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
    
    #==================================================
    #July
    #==================================================
    jul= ct.chart.line( month_cfs[6],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'Jul','yaxis': {'title':'CF(-)','range': [0,1],'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    jul= ct.chart.line( month_cfs_75[6],line=solar_q,fig= jul,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    jul= ct.chart.line( month_cfs_25[6],line=solar_q,fig=jul,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    jul= ct.chart.line( month_wfs[6],line=wind,fig= jul,scatter_kwargs={'mode':'lines','name':'wind' })   
    jul= ct.chart.line( month_wfs_75[6],line=wind_q,fig= jul,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    jul= ct.chart.line( month_wfs_25[6],line=wind_q,fig= jul,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
    
    
    #==================================================
    #August
    #==================================================
    aug= ct.chart.line( month_cfs[7],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'Aug','yaxis': {'title':'CF(-)','range': [0,1],'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    aug= ct.chart.line( month_cfs_75[7],line=solar_q,fig= aug,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    aug= ct.chart.line( month_cfs_25[7],line=solar_q,fig=aug,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    aug= ct.chart.line( month_wfs[7],line=wind,fig= aug,scatter_kwargs={'mode':'lines','name':'wind' })   
    aug= ct.chart.line( month_wfs_75[7],line=wind_q,fig= aug,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    aug= ct.chart.line( month_wfs_25[7],line=wind_q,fig= aug,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
    
    #==================================================
    #September
    #==================================================
    sep= ct.chart.line( month_cfs[8],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'Sep','yaxis': {'title':'CF(-)','range': [0,1],'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    sep= ct.chart.line( month_cfs_75[8],line=solar_q,fig= sep,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    sep= ct.chart.line( month_cfs_25[8],line=solar_q,fig=sep,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    sep= ct.chart.line( month_wfs[8],line=wind,fig= sep,scatter_kwargs={'mode':'lines','name':'wind' })   
    sep= ct.chart.line( month_wfs_75[8],line=wind_q,fig= sep,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    sep= ct.chart.line( month_wfs_25[8],line=wind_q,fig= sep,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
    
    #==================================================
    #October
    #==================================================
    octo= ct.chart.line( month_cfs[9],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'Oct','yaxis': {'title':'CF(-)','range': [0,1],'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    octo= ct.chart.line( month_cfs_75[9],line=solar_q,fig= octo,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    octo= ct.chart.line( month_cfs_25[9],line=solar_q,fig=octo,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    octo= ct.chart.line( month_wfs[9],line=wind,fig= octo,scatter_kwargs={'mode':'lines','name':'wind' })   
    octo= ct.chart.line( month_wfs_75[9],line=wind_q,fig= octo,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    octo= ct.chart.line( month_wfs_25[9],line=wind_q,fig= octo,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
    
    #==================================================
    #November
    #==================================================
    nov= ct.chart.line( month_cfs[10],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'Nov','yaxis': {'title':'CF(-)','range': [0,1],'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    nov= ct.chart.line( month_cfs_75[10],line=solar_q,fig= nov,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    nov= ct.chart.line( month_cfs_25[10],line=solar_q,fig=nov,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    nov= ct.chart.line( month_wfs[10],line=wind,fig= nov,scatter_kwargs={'mode':'lines','name':'wind' })   
    nov= ct.chart.line( month_wfs_75[10],line=wind_q,fig= nov,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    nov= ct.chart.line( month_wfs_25[10],line=wind_q,fig= nov,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
    
    #==================================================
    #December
    #==================================================
    dec= ct.chart.line( month_cfs[11],line=solar,scatter_kwargs={
       'mode':'lines','name':'solar'},layout_kwargs={'title':'Dec','yaxis': {'title':'CF(-)','range': [0,1], 'ticks':'outside','showline':True,
'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18},'xaxis':xaxis}) 
    dec= ct.chart.line( month_cfs_75[11],line=solar_q,fig= dec,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    dec= ct.chart.line( month_cfs_25[11],line=solar_q,fig=dec,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
    dec= ct.chart.line( month_wfs[11],line=wind,fig= dec,scatter_kwargs={'mode':'lines','name':'wind' })   
    dec= ct.chart.line( month_wfs_75[11],line=wind_q,fig= dec,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
    dec= ct.chart.line( month_wfs_25[11],line=wind_q,fig= dec,scatter_kwargs={'mode':'lines','name':'Q1'},fill='tonexty',
    fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
    
    

    #====================================================================
    #         Pearson Correlation Coeffieicnt 
    # resampling  hourly CF to  monthly data  and calulate peasron coeff.
    #==================================================================== 
    pearson=[]
    for i in range(len(year)):
            print('Caluclating Pearson cor. coeff.')
            print(year[i])
            CF_solar_cor=ct.cube.select(SCF_full, time=str(year[i]))
            CF_wind_cor=ct.cube.select(WCF_full, time=str(year[i]))
            #===================================================
            CFs=ct.cube.resample(CF_solar_cor,how='mean', freq='month')
            CFw=ct.cube.resample(CF_wind_cor, how='mean', freq='month')
            #===========================================+
            #calculating pearson cor 
            #============================================
            yearlymean= ct.cube.average(CFs, dim='time')
            sub_CFs= ct.operator.sub(CFs,yearlymean)
            sub_CFs= ct.cdm.update_attributes(sub_CFs,attrs={'units':'1'})
            yearlymean= ct.cube.average(CFw, dim='time')
            sub_CFw= ct.operator.sub(CFw,yearlymean)
            sub_CFw= ct.cdm.update_attributes(sub_CFw,attrs={'units':''})
            Num=ct.operator.mul(sub_CFs,sub_CFw)
            #print('Bug Area Passed')
            yearlysum=ct.cube.resample (Num, freq='year', dim='time', how='sum')
            CFs_sq=ct.math.square(sub_CFs)
            CFs_sq_sum=ct.cube.resample (CFs_sq, freq='year', dim='time', how='sum')
            CFw_sq=ct.math.square(sub_CFw)
            CFw_sq_sum=ct.cube.resample (CFw_sq, freq='year', dim='time', how='sum')
            #updating attributes to prevent bugs
            CFs_sq_sum= ct.cdm.update_attributes(CFs_sq_sum,attrs={'units':'1'})
            CFw_sq_sum= ct.cdm.update_attributes(CFw_sq_sum,attrs={'units':'1'})
            den_mul=ct.operator.mul(CFs_sq_sum,CFw_sq_sum)
            #print('Bug2 Area Passed')
            sqrt_den_mul=ct.math.sqrt(den_mul)
            pearson_r=ct.operator.div( yearlysum,sqrt_den_mul)
            #making pearson go from 0 to 1
            pearson_r=ct.operator.sub(1,pearson_r)
            pearson_r=ct.operator.mul(0.5,pearson_r)
            #print('Bug3 Area Passed')
            pearson.append(pearson_r)            
            
    #concatenating all pearson coeff. to full years
    for i in range(len(pearson)-1):
        pearson[i+1]=ct.cube.concat([pearson[i],pearson[i+1]],dim='time')
    pearson_full=pearson[-1]
    
    #taking avreage over time dimension  
    pearson_cor=ct.cube.average(pearson_full,dim='time')
    #pearson_cor=ct.cdm.update_attributes(pearson_cor,attrs={'units':''})
    pearson_cor=ct.cdm.update_attributes(pearson_cor, attrs={'long_name': 'r'})
    pearson_cor=ct.cube.select( pearson_cor, extent=[-180, 180, -60, 60])
    
    #==================================================================   
    # Plotting monthly solar and wind capacity factor profile
    #==================================================================   
    SCF_full_month=ct.cube.resample(SCF_full, how='mean', freq='month')
    WCF_full_month=ct.cube.resample(WCF_full, how='mean', freq='month')
        
    # x-axis properties 
    xaxis={'title':'Month','tickvals':[1,2,3,4,5,6,7,8,9,10,11,12],'ticktext':['J','F','M','A','M','J','J','A','S','O','N','D'],'ticks':'outside','showline':True,'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18}
    
    # checking year length 
    print('Checking year length:plot mean if year ==1')
    num_of_years=len(year)
    if num_of_years==1:
    
        #===========================================
        # Aggregating and selcting monthly solar CF
        #==========================================
        selcfs=ct.cube.groupby_reduce(SCF_full_month, group='time.month', how='mean', dim='time')
        selcfs=ct.geo.extract_point(selcfs, lon=lon, lat=lat)
        #======================================
        # Aggregating and selcting monthly wind
        #========================================
        selwfs=ct.cube.groupby_reduce(WCF_full_month, group='time.month', how='mean', dim='time')
        selwfs=ct.geo.extract_point(selwfs, lon=lon, lat=lat)
        #=====================
        # plots are made here
        #=======================
        mon_prof= ct.chart.line( selcfs,line=solar,scatter_kwargs={
            'mode':'lines','name':'CF_solar_mon'},layout_kwargs={'xaxis': xaxis,'yaxis': {'title':'CF(-)','range':[0,1],                          'ticks':'outside','showline':True,'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18}})
        mon_prof= ct.chart.line( selwfs,line=wind,fig=mon_prof,scatter_kwargs={'mode':'lines','name':'CF_wind_mon'})
        
    else: 
        # Quantiles are plottted here for year greater than 1
        #====================================================
        # Aggregating and selcting monthly solar CF
        #====================================================
        selcfs_med=ct.cube.groupby_reduce(SCF_full_month, group='time.month', how='quantile',q=0.5, dim='time')
        selcfs_q1=ct.cube.groupby_reduce(SCF_full_month, group='time.month', how='quantile',q=0.25, dim='time')
        selcfs_q3=ct.cube.groupby_reduce(SCF_full_month, group='time.month', how='quantile',q=0.75, dim='time')
        selcfs_med=ct.geo.extract_point(selcfs_med, lon=lon, lat=lat)
        selcfs_q1=ct.geo.extract_point(selcfs_q1, lon=lon, lat=lat)
        selcfs_q3=ct.geo.extract_point(selcfs_q3, lon=lon, lat=lat)
        
        #======================================
        # Aggregating and selcting monthly wind
        #=======================================
        selwfs_med=ct.cube.groupby_reduce(WCF_full_month, group='time.month', how='quantile',q=0.5, dim='time')
        selwfs_q1=ct.cube.groupby_reduce(WCF_full_month, group='time.month', how='quantile',q=0.25, dim='time')
        selwfs_q3=ct.cube.groupby_reduce(WCF_full_month, group='time.month', how='quantile',q=0.75, dim='time')
        selwfs_med=ct.geo.extract_point(selwfs_med, lon=lon, lat=lat)
        selwfs_q1=ct.geo.extract_point(selwfs_q1, lon=lon, lat=lat)
        selwfs_q3=ct.geo.extract_point(selwfs_q3, lon=lon, lat=lat)
        #=========================================
        # plots and qauntiles are made here
        #=========================================
        mon_prof= ct.chart.line( selcfs_med,line=solar,scatter_kwargs={'mode':'lines','name':'CF_solar_mon'},layout_kwargs={'xaxis': xaxis,'yaxis': {'title':'CF(-)','range':[0,0.6],      'ticks':'outside','showline':True,'linewidth':1,'linecolor':'black','mirror':True,'tickfont_size':18}})
        mon_prof= ct.chart.line( selcfs_q3,line=solar_q,fig= mon_prof,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
        mon_prof= ct.chart.line(selcfs_q1,line=solar_q,fig=mon_prof,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(255, 0, 0, 0.2)',showlegend=False)
        
        mon_prof= ct.chart.line( selwfs_med,line=wind,fig=mon_prof,scatter_kwargs={'mode':'lines','name':'CF_wind_mon'})
        mon_prof= ct.chart.line( selwfs_q3,line=wind_q,fig= mon_prof,scatter_kwargs={'mode':'lines','name':'Q3' },showlegend=False)
        mon_prof= ct.chart.line(selwfs_q1,line=wind_q,fig=mon_prof,scatter_kwargs={'mode':'lines','name':'Q1' },fill='tonexty',fillcolor='rgba(0, 0, 255, 0.2)',showlegend=False)
        

   #==================================================================
    #    Annual average cstab & pearson coeff. 
    # Please note that annual cstab & pearson coeff. are printed to   
    # the console which are then used to calculate the standard deviation
    # of capicity factors offline 
    #==================================================================
    print('Annual average,std cstab & pearson')
    #Mean Cstab
    an_avg_cstab=ct.geo.extract_point(Cstab_Avg, lon=lon, lat=lat)
    an_avg_cstab=ct.cdm.get_value(an_avg_cstab)
    print('Annual average cstab='+str(an_avg_cstab))
    #Mean pearson 
    an_avg_pc=ct.geo.extract_point( pearson_cor, lon=lon, lat=lat)
    an_avg_pc=ct.cdm.get_value(an_avg_pc)
    print('Annual average pearson='+str(an_avg_pc))
    
    #standard deviation cstab
    Cstab_an_full=[]
    for i in range(len(year)):
        Cstab_an=ct.cube.select(Cstab_full, time=str(year[i]))
        Cstab_an_av=ct.cube.average(Cstab_an, dim= 'time')
        cstab_forstd=ct.geo.extract_point(Cstab_an_av, lon=lon, lat=lat)
        Cstab_an_full.append(cstab_forstd)
    print('Annual average data for STD of cstab is')
    print(Cstab_an_full)
      
    #standard deviation pearson
    pearson_an_full=[]
    for i in range(len(year)):
        pearson_an=ct.cube.select(pearson_full, time=str(year[i]))
        pearson_an_av=ct.cube.average(pearson_an, dim='time')
        pearson_forstd=ct.geo.extract_point(pearson_an_av, lon=lon, lat=lat)
        pearson_an_full.append(pearson_forstd)
    print('Annual average data for STD of pearson corr. is') 
    print(pearson_an_full)
      
    #==================================================================
    # Maps_monthly
    #==================================================================
    #update pearson title and attributes 
    MAP_CONFIG['title']['text_lines']=["Seasonal Synergy (-)"]
    MAP_CONFIG['contour'][ 'contour_interpolation_ceiling']=ct.math.nanmax(pearson_cor)
    MAP_CONFIG['contour']['contour_max_level']=ct.math.nanmax(pearson_cor)
   
    fig_Pearson = ct.map.plot(pearson_cor,**MAP_CONFIG)
    
    print('max_pearson_cor_'+str(ct.math.nanmax(pearson_cor)))
    print('max_Cstab_Avg_'+str(ct.math.nanmax(Cstab_Avg)))
    
    #=============================================
    #Extrating threshold  for good synergy as float
    #==============================================
    T_hold=[]
    for i in thres.split('/'):
        T_hold.append(float(i))
    print(T_hold)
    pearson_cor_good=ct.cube.where(ct.math.logical_not                 (ct.math.logical_and(pearson_cor>T_hold[2],Cstab_Avg>T_hold[0])),pearson_cor,9)#good daily, good seasonal
    
    pearson_cor_good=ct.cube.where(ct.math.logical_not                 (ct.math.logical_and(ct.math.logical_and(pearson_cor>T_hold[3],pearson_cor<T_hold[2]),Cstab_Avg>T_hold[0])),pearson_cor_good,8)#good daily, medium seasonal
         
    pearson_cor_good=ct.cube.where(ct.math.logical_not                 (ct.math.logical_and(pearson_cor<T_hold[3],Cstab_Avg>T_hold[0])),pearson_cor_good,7)#  good daily, bad seasonal
    
    pearson_cor_good=ct.cube.where(ct.math.logical_not                 (ct.math.logical_and(pearson_cor>T_hold[2],ct.math.logical_and(Cstab_Avg>T_hold[1],Cstab_Avg<T_hold[0]))),pearson_cor_good,6)# medium daily, good seasonal
    
    pearson_cor_good=ct.cube.where(ct.math.logical_not                 (ct.math.logical_and(ct.math.logical_and(pearson_cor>T_hold[3],pearson_cor<T_hold[2]),ct.math.logical_and(Cstab_Avg>T_hold[1],Cstab_Avg<T_hold[0]))),pearson_cor_good,5)# medium daily, medium seasonal
    
    pearson_cor_good=ct.cube.where(ct.math.logical_not                 (ct.math.logical_and(pearson_cor<T_hold[3],ct.math.logical_and(Cstab_Avg>T_hold[1],Cstab_Avg<T_hold[0]))),pearson_cor_good,4)# medium daily, bad seasonal
    
    pearson_cor_good=ct.cube.where(ct.math.logical_not                 (ct.math.logical_and(pearson_cor>T_hold[2],Cstab_Avg<T_hold[1])),pearson_cor_good,3)# bad daily, good seasonal 
    
    pearson_cor_good=ct.cube.where(ct.math.logical_not                 (ct.math.logical_and(ct.math.logical_and(pearson_cor>T_hold[3],pearson_cor<T_hold[2]),Cstab_Avg<T_hold[1])),pearson_cor_good,2)# bad daily, medium seasonal 
    
    pearson_cor_good=ct.cube.where(ct.math.logical_not                 (ct.math.logical_and(pearson_cor<T_hold[3],Cstab_Avg<T_hold[1])),pearson_cor_good,1)# bad daily, bad seasonal 
    
    
    
 
    
    pearson_cor_good=ct.cube.where(SCF>=0.1,pearson_cor_good)# low solar resource 
    pearson_cor_good=ct.cube.where(WCF>=0.1,pearson_cor_good)# low wind resource
    
    
    #=================================================================
    # plotting Good and bad synergy zones
    #==================================================================
    MAP_CONFIG['title']['text_lines']=["Synergy Zones"]
    MAP_CONFIG['contour']['contour_shade_colour_method']="list"
    MAP_CONFIG['contour']['contour_shade_colour_list'] =[
"rgb(76,0,83)","rgb(75,44,123)","rgb(48,83,139,)", "rgb(0,116,143)","rgb(0,147,142)", "rgb(0,177,133)", "rgb(0,204,106)", "rgb(159,221,63)", "rgb(255,230,42)",]
    MAP_CONFIG['contour']['contour_interpolation_ceiling']=10
    MAP_CONFIG['contour']['contour_interpolation_floor']="off"
    MAP_CONFIG['contour']['contour_min_level']="off"
    MAP_CONFIG['contour']['contour_max_level']="off"

    fig_Pearson_good= ct.map.plot(pearson_cor_good,**MAP_CONFIG)
         
    return fig_WCF,fig_SCF,fig_Cstab,fig_Pearson, fig_Pearson_good,Cstab_ma,mon_prof,jan,feb,mar,apr,may,jun,jul,aug,sep,octo,nov,dec
