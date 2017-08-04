function CCaPS
%
clear
clc
% Box model is based on the models of Kump and Arthur 1999 and Payne and
% Kump 2007. Carbonate system solver is based on Emerson and Hedges 2008
% implementation of R. Zeebe's Csys code (can be found on his website). 
%
% Model sketch:
%                                   
%                       _______
%          Fwcarb-->   |       |--> Fbcarb 
%          Fworg -->   |  MC   |
%          Fvolc  -->  |       |--> Fborg 
%                      |_______|
%                       _______
%         Fwcarb-->    |       |
%         Fwsil -->    | MCa   |   --> Fbcarb
%         Fhyd  -->    |       |                                
%                      |_______|
% 
%                       _______ 
%         Fwp   -->    |  MP   |   --> Fbp                                
%                      |_______|   
%
%                       _______
%         Fborg -->    |  O2   |   --> Fworg + Fvolc_reduced                               
%                      |_______|   
% 
% 
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%
% Choose run type HERE 
% Uncomment to re-create paper figures or do a 'custom' run. 
% If the latter set your own preferences in custom_params function

% flags.run = 'LBC_plot';     %%(Figure 1)
% flags.run = 'short_pert';   %%(Figure 3)
% flags.run = 'long_pert';    %%(Figure 4)
%  flags.run = 'Mp_sens';      %%(Figure 5)
%flags.run = 'multi_pert';   %%(Figure 6)
flags.run = 'custom';

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch flags.run
    case 'short_pert'
        Output = short_pert(flags);
    case 'long_pert'
       Output = long_pert(flags);
   case 'Mp_sens'
       Output = Mp_sens(flags);
    case 'multi_pert'
        Output = multi_pert(flags);
    case 'custom'      
        Output = custom_params(flags);
    case 'LBC_plot'      
        LBC_plot();
        Output = [];
    otherwise
        error('No such run option!')
end

%export model vars to base workspace
assignin('base','Output',Output)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Output] = custom_params(flags)
   

   

%run length [yrs] 
params.t_run_start = -100e3;
params.t_run_end = 2500e3;



%flag0 - pco2 basline
flags.pCO2 = 4;    
                    %0 1.3 ppm (alk=.9, dic=.3)
                    %1 502 ppm (alk=2.2, dic=2.0)
                    %2 37000 ppm (alk=19.0, dic=20.0)
                    %3 2000 ppm (alk=5.0, dic=5.22 omega = 9)
                    %4 2000 ppm (alk =1.89, dic = 1.9 omega  = 1.4)

%flag1 - pco2 weathering dependence:
flags.weathering = 2;     
                    %0 - off
                    %1 - on, dependent on RCO2 with a power of 
                    params.p_w = 1; %1 is linear, 2 is cubic, etc..
                    %2 - on, Berner parameterization



%flags.del_diff - differences between delccarb and delcorg:
flags.del_diff = 3;      
                    %0 - off
                    %1 - extrabasinal Corg. 
                        %amount: 
                        params.fexternal = .5; 
                        %composition:
                        params.delexternal = -28.0;
                    %2 - pco2 depd. fractionation with const PO4
                    %3-pco2 depd. fractionation with variable PO4 using 
                       %parameterization from Kump and Arthur 99
                       %note that might not yield baseline of +2
                       


%flags.Fbcarb - Fbcarb feedback                      
flags.Fbcarb = 5;         
                    %0- Fbcarb = const initial Fbcarb
                    %1- Fbcarb= inputs of the Ca box 
                    %2- Fbcarb = inputs of C box minus Fborg.
                    %3- linearly proportional to [Ca]
                    %4- proportional to omega with a power of 
                        params.p_carb = 1;    
                    %5- proportional to omega with negative omega<1    


                           


%flags.pertb - perturbtion type                    
flags.pertb = 3.4;       
                    
    %0 - steady state

    %1 - drive a gaussian change in Fborg
    params.f_extra = .66; %additional multiples of initial 
    %Fborg. 


    %3.x (see below) change in the carbon input (Fextra); 
    %perturbation length - must be shorter than run length!!
    %and no overlap between perturbations (flag 3.9)
    
    
    params.t_perturbation_start_1 = 0;
    params.t_perturbation_end_1 = 50e3;     
    
    params.M_extra_1= 3.5e17; %total additional carbon in moles
%     params.M_extra_1= 1e18;
%     params.del_extra_1 = -70;  %isotopic composition
    params.del_extra_1 = -5.5; 
    

    %3.0 %Fborg is the C inputs minus Fbcarb (not to be
         %used with flags.Fbcarb = 1 )
                                    %
    %3.1 %Fborg is linearly proportional to RCO2

    %3.2 %Fborg is prop to P burial; C:P is set at 106 (Redfield)

    %3.3 %Fborg is prop to P burial;C:P varies as 
         %a func of pco2 with a power of
         params.p_CP_RCO2=1;

    %3.4 %Fborg remains at the original steady
    %state

    %3.5 %Fborg is such that forg remains
    %unchanged

    %3.6 %Fborg is prop to P burial; C:P varies
    %as a power of Fbp with a power of 
     params.p_CP_bP = 1;

    %3.7 %%Fborg is prop to P burial; C:P varies
    %as an exponential func of Fbp


    %3.8 Fborg is prop to p weathering, C:P at
    %106

    %3.9 Add a second perturbation
    % No overlap between the two perturbations!!! 
    % (Fborg as in 3.6).

    params.t_perturbation_start_2 = 20e3;

    params.t_perturbation_end_2 = 800e3;     

    params.M_extra_2 = 1.2e18;

    
    params.del_extra_2 = -5.5;


    %4 Add another 2 pulses (4 total)

    params.t_perturbation_start_3 = 270e3 + 100e3;

    params.t_perturbation_end_3 = 300e3 + 270e3 + 100e3;     

    params.M_extra_3 = 0e18;

    params.del_extra_3 = -5.5;

    

    params.t_perturbation_start_4 = 620e3 + 100e3;

    params.t_perturbation_end_4 = 100e3 + 620e3 + 100e3;     

    params.M_extra_4 = 0e17;

    params.del_extra_4 = -5.5;
                            
    
    %5 add another pulse and a shift in silicate weathering
    
    params.t_perturbation_start_5 = 840e3 + 100e3;
    
    params.t_perturbation_end_5 = 100e3 + 840e3 + 100e3;     

    params.M_extra_5 = 7.5e17;
    
    params.del_extra_5 = -5.5;
    
    
    params.t_shift_start = 500e3 + 100e3;

    params.shift_mag = 3e12;
    
    
    %case 5.1- shift to a different C:P power
    params.p_CP_bP_2 = 1;
    
    params.t_p_CP_bP = 600e3;

%flags.PO4 - phophate burial
flags.PO4 = 0;          %0 - phosphate burial is prop to phophate mass
                    %with a power of 
                    params.p_bp_Mp = 1; 

params.PO4_multi = 1; %multiplyer of initial PO4 (set at modern)
                    

%flags.plot - number and type of plotted results. 
flags.plot = 6; 
                    %0- plot all 20 output parameters: 
                        %1 DIC;
                        %2 d13C carb;
                        %3 Ca;
                        %4 Fborg
                        %5 Fbcarb;
                        %6 forg;
                        %7 pCO_2;
                        %8 pH;
                        %9 Omega;
                        %10 CO2aq;
                        %11 HCO3;
                        %12 avg_d13C of weathering;
                        %13 ALK; 
                        %14 average d13C of burial;
                        %15 d13C org';
                        %16 PO4;
                        %17 d44Ca;
                        %18 Fextra;
                        %19 C:P;
                        %20 Fwsil;
                        %21 CO2
                
                    %1- plot only 9 outputs: 
                       %Fextra, d13C, Ca, DIC, pCO2, ALk, d44Ca, C:P, Omega
                       
                    %2- plot only 6 outputs:
                       %Fextra, d13C, pCO2, forg, Omega, C:P
                      
                    %2.1- plot only 6 outputs:
                       %Fextra, d13C, pCO2, forg, Omega, C:P  
                    
                    %2.2- plot only 6 outputs:
                       %Fextra, d13C, pCO2, forg, Omega, pH     
                       
                    %3- plot only 4 outputs:
                        %forg, d13C, ALK, pCO2
                    
                    %3.1- plot only 4 outputs:
                        %Fextra, d13C, Omega, pco2
                    
                    %4- plot only 3 outputs:
                        %Fborg, d13C, Omega
                        
                        
                    %5- plot only 2 outputs:
                        %Fextra, d13C
                    
                    %5.1 plot only 2 outputs:
                        %d13C, Omega
                        
                        
                    %6- plot only delorg and delcarb on the same plot 
                    
                    %7- plot trajectory on alk-dic phase plot
                    

%flags.plot_save - plot save options    
flags.plot_save = 0;
                    %0- don't save a plot
                    
                    %1- save a copy as an eps file. If putting file in
                    %directory other than the current one then add the full
                    %path to the filename. 
                    

      flags.file_name  = 'myfile';

      flags.file_path  = '/some/very/long/path';               
                        
      

%flags.plot_hold - plot overlay      
flags.plot_hold = 'hold';      %hold - do not clear plot window
                        %replace - clear window
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                       
    
    %call the model
    [Output] = model(params,flags);
 
    %call plotting function
    plotting(Output,params,flags);
    
 
end

function [] = LBC_plot()

close all

hf = figure('Units', 'centimeters','OuterPosition', ...
           [1 1 19 16],'name','sim_window','PaperPositionMode','auto');  

[LBC_Data] = LBC_data();


ax1 = subplot('Position',[0.18 0.1 0.3 0.8],'FontSize',9);

LBC_carb_lines(LBC_Data)

ylabel('\delta^{13}C_{carb}');

ylim([-5 7])

set(gca,'YDir','reverse')

grid on

box on

view(-90,90) % swap the x and y axis

axes('Position',get(ax1,'Position'),...
'Visible','on','Color','none','XTickLabel',[],'XTick',[],'XColor','k',...
'YMinorTick','on','YAxisLocation','right','YDir','reverse','FontSize',9);

LBC_org_lines(LBC_Data)

ylabel('\delta^{13}C_{org}','VerticalAlignment','Middle',...
    'HorizontalAlignment','center','FontSize',9);
    
ylim([-35 -23])

view(-90,90) % swap the x and y axis

%%%%%%

ax2 = subplot('Position',[0.69 0.1 0.3 0.8],'FontSize',9);

LBC_carb_lines(LBC_Data)

ylabel('\delta^{13}C_{carb}','FontSize',9);

xlim([25 100])
ylim([-5 7])

grid on

box on

set(gca,'YDir','reverse')

view(-90,90) 

axes('Position',get(ax2,'Position'),...
'Visible','on','Color','none','XTickLabel',[],'XTick',[],'XColor','k',...
'YMinorTick','on','YAxisLocation','right','YDir','reverse','FontSize',9);


LBC_org_lines(LBC_Data)

ylabel('\delta^{13}C_{org}','VerticalAlignment','Middle',...
    'HorizontalAlignment','center','FontSize',9);

xlim([25 100]) 
ylim([-35 -23])

view(-90,90) 

% print('-f1','-depsc','-loose', '-painters', '-noui',...
% ['/Users/avivbachan/Documents/Research/Thesis/Modeling_chapter/Modelling_paper/Manuscript_version_2/Figures'...
% '/LBC_mat_fig.eps']);


end

function [Output] = short_pert(flags)   


params = struct(...
't_run_start',  -10e3,...
't_run_end',    200e3,...
'p_carb',       1,... 
'p_CP_bP',      1,...
'p_bp_Mp',      1,...
'PO4_multi',    1,...
't_perturbation_start_1', 0,...
't_perturbation_end_1',  20e3 ,...
'M_extra_1', 3.5e17,...
'del_extra_1', {-5.5, -22, -50, -70 }...
);



flags.pCO2 = 4;

flags.weathering = 2;

flags.del_diff = 0;

flags.Fbcarb = 5;

flags.pertb = 3.6;


flags.PO4 = 0;


flags.plot = 3.1;

flags.plot_save = 0; 

flags.plot_hold = 'hold';

flags.file_name  = '';

flags.file_path  = '';     
%%%%%%%%%%%%%%%%%



for i = 1:length(params)
    
    %call the model
    
    Output = model(params(i),flags);
 
    
    %call plotting function
    [hf,hs] = plotting(Output,params,flags);
end

    set(hs, 'XLim',[params(1).t_run_start, params(1).t_run_end ])

    %line colors black
    hline = findobj(hf, 'type', 'line');
    set(hline,'color','k')

    %dashed lines in the isotope panel
    h_fig_children = get(hf,'Children');
    h_delta_children = get(h_fig_children(3),'Children');
    h_delta_children = flipud(h_delta_children);
    set(h_delta_children(1),'LineStyle','-')
    set(h_delta_children(2),'LineStyle','--')
    set(h_delta_children(3),'LineStyle','-.')
    set(h_delta_children(4),'LineStyle',':')
    
    
    
    
    %legend
    legend(h_delta_children,...
        num2str([params.del_extra_1]'),'Location','SouthEast')
    
    %letters
    % Create textbox
    annotation(hf,'textbox',[0.06556 0.9027 0.06261 0.06701],...
    'String',{'(A)'},...
    'LineStyle','none');

    % Create textbox
    annotation(hf,'textbox',[0.5235 0.9027 0.06261 0.06701],'String',...
    {'(B)'},'LineStyle','none');

    % Create textbox
    annotation(hf,'textbox',[0.5235 0.4622 0.06261 0.06701],'String',...
    {'(D)'},'LineStyle','none');

    % Create textbox
    annotation(hf,'textbox',[0.06556 0.4622 0.06261 0.06701],...
    'String',{'(C)'},'LineStyle','none');
    
    



end

function [Output] = long_pert(flags)

params = struct(...
't_run_start',  -100e3,...
't_run_end',    800e3,...
'p_carb',       1,... 
'p_CP_bP',      1,...
'p_bp_Mp',      1,...
'PO4_multi',    1,...
't_perturbation_start_1', 0 ,...
't_perturbation_end_1',  {100e3 , 500e3  },...
'M_extra_1', {1.5e18  3e18},...
'del_extra_1', -5.5,...
'p_w' , 1,...
'p_CP_RCO2', 1);

params(2,:) = deal(params(1,:));
[params(2,:).p_CP_bP] = deal(3);

flags.pCO2 = 4;
flags.weathering = 2;
flags.del_diff = 0;
flags.Fbcarb = 5;
flags.pertb = 3.6;
flags.PO4 = 0;
flags.plot = 3.1;
% flags.plot = 0;
flags.plot_save = 0;     
flags.plot_hold = 'hold';
flags.file_name  = '';
flags.file_path  = '';               





%%%%%%%%%%%%%%%%

[m,n] = size(params);

for j = 1:n

for i = 1:m
    
    %call the model
    
    [Output] = model(params(i,j),flags);
 
    %call plotting function
    [hf,hs] = plotting(Output,params(i,j),flags);
end

end
    
%plot annotation

set(hs, 'XLim',[params(1).t_run_start, params(1).t_run_end ])

h_fig_children = get(hf,'Children');
h_lines_cell = get(h_fig_children,'Children');
h_lines_mat = [h_lines_cell{:}];
set(h_lines_mat,'Color','k') %black lines
set(h_lines_mat([1 3],1:3),'LineStyle','--','Color','k')
% % set(h_lines_mat(3,:),'LineStyle','-','Color','k')
% 
legend(h_lines_cell{1},'3','1',...
    'Location','NorthEast')
% 
% 
annotation(hf,'textbox',[0.06556 0.9027 0.06261 0.06701],...
'String',{'(A)'},...
'LineStyle','none');

% Create textbox
annotation(hf,'textbox',[0.5235 0.9027 0.06261 0.06701],'String',...
{'(B)'},'LineStyle','none');

% Create textbox
annotation(hf,'textbox',[0.5235 0.4622 0.06261 0.06701],'String',...
{'(D)'},'LineStyle','none');

% Create textbox
annotation(hf,'textbox',[0.06556 0.4622 0.06261 0.06701],...
'String',{'(C)'},'LineStyle','none');
%     

end

function [Output] = Mp_sens(flags)

params = struct(...
't_run_start',  - 100e3,...
't_run_end',    2000e3,...
'p_carb',       1,... 
'p_bp_Mp',      1,...
't_perturbation_start_1', 0 ,...
't_perturbation_end_1',  500e3,...
'M_extra_1', 2e18,...
'del_extra_1', -5.5,...
'p_w' , 1,...
'p_CP_RCO2', 1,...
'p_CP_bP', 1 ,...
'PO4_multi', 1);

params(2,:) = deal(params(1,:));
[params(2,:).PO4_multi] = deal(10);

params(3,:) = deal(params(2,:));
[params(3,:).PO4_multi] = deal(50);



flags.pCO2 = 4;
flags.weathering = 2;
flags.del_diff = 0;
flags.Fbcarb = 5;
flags.pertb = 3.6;
flags.PO4 = 0;
flags.plot = 3.1;
% flags.plot = 0;
flags.plot_save = 0;     
flags.plot_hold = 'hold';
flags.file_name  = '';
flags.file_path  = '';    

           




%% %%%%%%%%%%%%%%%%

[m,n] = size(params);

for j = 1:n

for i = 1:m

    
    %call the model
    
    [Output] = model(params(i,j),flags);
 
    %call plotting function
    [hf,hs] = plotting(Output,params(i,j),flags);

    
end

end


set(hs, 'XLim',[params(1).t_run_start, params(1).t_run_end ])

    
%plot annotation

h_fig_children = get(hf,'Children');
h_lines_cell = get(h_fig_children,'Children');
h_lines_mat = flipud([h_lines_cell{:}]);
set(h_lines_mat(1,:),'LineStyle','-.','Color','k')
set(h_lines_mat(2,:),'LineStyle','--','Color','k','LineWidth',1)
set(h_lines_mat(3,:),'LineStyle','-','Color','k')


%legend
legend(h_lines_mat(:,1),...
    '1x',...
    '10x', ...
    '50x','Location','NorthEast')




annotation(hf,'textbox',[0.06556 0.9027 0.06261 0.06701],...
'String',{'(A)'},...
'LineStyle','none');

% Create textbox
annotation(hf,'textbox',[0.5235 0.9027 0.06261 0.06701],'String',...
{'(B)'},'LineStyle','none');

% Create textbox
annotation(hf,'textbox',[0.5235 0.4622 0.06261 0.06701],'String',...
{'(D)'},'LineStyle','none');

% Create textbox
annotation(hf,'textbox',[0.06556 0.4622 0.06261 0.06701],...
'String',{'(C)'},'LineStyle','none');
    

end

function [Output] = multi_pert(flags)

params.t_run_start = -100e3;
params.t_run_end = 1000e3;

flags.pCO2 = 4;
flags.weathering = 2;
flags.del_diff = 0;
flags.Fbcarb = 5;

flags.pertb = 5.1;


flags.PO4 = 0;

% flags.plot = 5.3; %model d13C and Lombardy Basin d13C (no corr lines)
% flags.plot = 5.4; %Fextra, model d13C and Lombardy basin d13C
% flags.plot = 5.5; %pCO2, model d13C and Lombardy basin d13C
% flags.plot = 5.2; %Omega and PCO2
% flags.plot = 3.2; %C:P, d13C, Omega, pCO2 
flags.plot = 3.1; %Fextra, d13C, Omega, pCO2 
% flags.plot = 1;

flags.plot_save = 0;     
flags.plot_hold = 'hold';

%params used
params.p_carb = 1;
params.p_CP_bP = 3;
params.p_bp_Mp = 1;
params.PO4_multi = 1;




params.t_perturbation_start_1 = 0 ;
params.t_perturbation_end_1 = 20e3 ;
params.M_extra_1= 3.5e17; %total additional carbon in moles
params.del_extra_1 = -70;

params.t_perturbation_start_2 = 10e3 ;
params.t_perturbation_end_2 = 70e3 ;     
params.M_extra_2 = 1.2e18;
params.del_extra_2 = -5.5;

params.t_perturbation_start_3 = 270e3 ;
params.t_perturbation_end_3 = 150e3 + 270e3 ;     
params.M_extra_3 = 2e18;
params.del_extra_3 = -5.5;


params.t_perturbation_start_4 = 620e3 ;
params.t_perturbation_end_4 = 720e3 ;     
params.M_extra_4 = 7.5e17;
params.del_extra_4 = -5.5;

params.t_perturbation_start_5 = 840e3 ;
params.t_perturbation_end_5 = 940e3 ;     
params.M_extra_5 = 7.5e17;
params.del_extra_5 = -5.5;

params.t_shift_start = 400e3 ;
params.shift_mag = 3e12;

params.p_CP_bP_2 = 1;
params.t_p_CP_bP = 500e3;



 
%not in use 
params.p_w = 1;
params.fexternal = .5; 
params.delexternal = -28.0;
params.f_extra = .66;
params.p_CP_RCO2=1;
flags.file_name  = 'myfile';
flags.file_path  = '/some/very/long/path';               




%%%%%%%%%%%%%%%%%

    
    %call the model
    
 [Output] = model(params,flags);
    
 [hf, hs] = plotting(Output,params,flags);
   
switch flags.plot
    
    case 5.2
         
    set(hs,'XLim',[params.t_run_start params.t_run_end])
     
    case {5.3,5.4,5.5}
        
        
        switch flags.plot
            case {5.4,5.5}

        set(hs(2),'Position',[0.13 0.48 0.775 0.3],...
            'YLim',[-5 7])
        
        switch flags.plot
            case 5.4
        set(hs(1), 'Position',[0.13 0.8 0.775 0.1],...
            'XTickLabel','')
            case 5.5
                set(hs(1), 'Position',[0.13 0.8 0.775 0.1],...
            'XTickLabel','','YLim',[1000 5000])
        end
        
        
       
        set(get(hs(1),'XLabel'),'String','')
        
        
        annotation(hf,'textbox',[0.04619 0.9169 0.06261 0.06701],...
        'String',{'(A)'},...
        'FitBoxToText','off',...
        'LineStyle','none');

        annotation(hf,'textbox',[0.04619 0.67351 0.06261 0.06701],...
        'String',{'(B)'},...
        'FitBoxToText','off',...
        'LineStyle','none');
    
    
        annotation(hf,'textbox',[0.04619 0.3607 0.06261 0.06701],...
        'String',{'(C)'},...
        'FitBoxToText','off',...
        'LineStyle','none');
    
    %%%%%
        % Create line
annotation(hf,'line',[0.1975 0.1975],[0.8968 0.09817],'LineStyle','--',...
    'LineWidth',1);

% Create line
annotation(hf,'line',[0.2228 0.2222],[0.395 0.09651],'LineStyle','--',...
    'LineWidth',1);

% Create line
annotation(hf,'line',[0.3668 0.2205],[0.4765 0.3977],'LineStyle','--',...
    'LineWidth',1);

% Create line
annotation(hf,'line',[0.3686 0.3686],[0.8952 0.4775],'LineStyle','--',...
    'LineWidth',1);

% Create line
annotation(hf,'line',[0.5038 0.5032],[0.3966 0.09814],'LineStyle','--',...
    'LineWidth',1);

% Create line
annotation(hf,'line',[0.5897 0.5044],[0.4782 0.3943],'LineStyle','--',...
    'LineWidth',1);

% Create line
annotation(hf,'line',[0.5914 0.5914],[0.8968 0.4792],'LineStyle','--',...
    'LineWidth',1);

% Create line
annotation(hf,'line',[0.6473 0.6467],[0.3982 0.09977],'LineStyle','--',...
    'LineWidth',1);

% Create line
annotation(hf,'line',[0.7331 0.6479],[0.4798 0.396],'LineStyle','--',...
    'LineWidth',1);

% Create line
annotation(hf,'line',[0.7349 0.7349],[0.8984 0.4808],'LineStyle','--',...
    'LineWidth',1);
    
    
        end
        
        
        switch flags.plot
            case 5.3
            set(gca,'Position',[0.13 0.6 0.77 0.35],'YLim',[-5 7]);
            
            annotation(hf,'textbox',[0.04619 0.9169 0.06261 0.06701],...
            'String',{'(A)'},...
            'FitBoxToText','off',...
            'LineStyle','none');

            % Create textbox
            annotation(hf,'textbox',[0.04619 0.4622 0.06261 0.06701],...
            'String',{'(B)'},...
            'FitBoxToText','off',...
            'LineStyle','none');
            
        end
        
        



if get(hf,'UserData')==1 %if this is the first time plotted

%get LBC data    
LBC_Data = LBC_data();

%carbonate carbon plotting

switch flags.plot
    case {5.4,5.5}
     ax1 = axes('Position',[0.13 0.1 0.77 0.3]);
end
      
      
switch flags.plot
    case 5.3
     ax1 = axes('Position',[0.13 0.1 0.77 0.35]);
end


LBC_carb_lines(LBC_Data)

ylabel('\delta^{13}C_{carb}');

ylim([-5 7])

grid on

box on


%organic carbon plotting

ax2 =  axes('Position',get(ax1,'Position'),...
   'Visible','on','Color','none','XTickLabel',[],'XTick',[],'XColor','k',...
    'YMinorTick','on','YAxisLocation','right');


LBC_org_lines(LBC_Data)


ylabel('\delta^{13}C_{org}');

xlabel('Thickness','Color','k','Position',[297.86 -36.949 17.321]);
    
ylim([-35 -23])

end
    
    case 3.1
   
        set(hs, 'XLim',[params(1).t_run_start, params(1).t_run_end ])
        
        h_fig_children = get(hf,'Children');
        h_lines_cell = get(h_fig_children,'Children');
        h_lines_mat = flipud([h_lines_cell{:}]);
        set(h_lines_mat(:),'Color','k')
        
        annotation(hf,'textbox',[0.06556 0.9027 0.06261 0.06701],...
        'String',{'(A)'},...
        'LineStyle','none');

        % Create textbox
        annotation(hf,'textbox',[0.5235 0.9027 0.06261 0.06701],'String',...
        {'(B)'},'LineStyle','none');

        % Create textbox
        annotation(hf,'textbox',[0.5235 0.4622 0.06261 0.06701],'String',...
        {'(D)'},'LineStyle','none');

        % Create textbox
        annotation(hf,'textbox',[0.06556 0.4622 0.06261 0.06701],...
        'String',{'(C)'},'LineStyle','none');
        
end


end

function [hf,hs]= plotting(Output,params,flags)
%plotting function. called from control function
    

%unpackage the model output
x=Output.t;
y=[Output.DIC, Output.delbcarb, Output.Ca ,...
   Output.Fborg, Output.Fbcarb,...
   Output.forg, Output.pco2, Output.pH,...
   Output.omega, Output.epsb, Output.hco3,...
   Output.delwC, Output.ALK, Output.delbC,...
   Output.delborg,Output.PO4,Output.delCa,...
   Output.Fextra,Output.CP,Output.Fwsil,...
   Output.CO2];




%open a new figure window if one isnt already opened
if isempty(findobj('type','figure','name','sim_window'))==1,
    close all
    
    switch flags.plot
        case 7
           hf = figure('name','sim_window'); 
        case 6
           fig_height = 19;
           hf = figure('Units', 'centimeters','OuterPosition', ...
               [1 1 19 fig_height],'name','sim_window'); 
        case {5,5.1}
           fig_height = 9;
           hf = figure('Units', 'centimeters','OuterPosition', ...
               [1 1 19 fig_height],'name','sim_window'); 
       case 5.2
           fig_height = 16;
           hf = figure('Units', 'centimeters','OuterPosition', ...
               [1 1 8 fig_height],'name','sim_window'); 
       case 5.3
           fig_height = 16;
           hf = figure('Units', 'centimeters','OuterPosition', ...
               [1 1 16 fig_height],'name','sim_window');     
       case {5.4,5.5}
           fig_height = 19;
           hf = figure('Units', 'centimeters','OuterPosition', ...
               [1 1 16 fig_height],'name','sim_window');         
           
       case 4
           fig_height = 7;
           hf = figure('Units', 'centimeters','OuterPosition', ...
               [1 1 19 fig_height],'name','sim_window'); 
       case {3,3.1,3.2,3.3}
           fig_height = 13;
           hf = figure('Units', 'centimeters','OuterPosition', ...
               [1 1 16 fig_height],'name','sim_window'); 
        case 2
            fig_height = 20;
            hf = figure('Units', 'centimeters','OuterPosition', ...
                [1 1 19 fig_height],'name','sim_window'); 
        case {2.1,2.2}
            fig_height = 12;
            hf = figure('Units', 'centimeters','OuterPosition', ...
                [1 1 19 fig_height],'name','sim_window'); 
        case {1,1.1,0}
            hf = figure('Units', 'normalized','OuterPosition', ...
                [0 0 1 1],'name','sim_window'); 
    end
    
    fig_counter = 1;
    
    set(hf,'UserData',fig_counter);
    
else 
    
    switch flags.plot_hold
        case 'replace'
        clf;
        fig_counter = 1;
        case 'hold'
            hf = gcf;
            set(hf,'NextPlot','add')
            fig_counter = get(gcf,'UserData') + 1;
            set(hf,'UserData',fig_counter);
    end
end
   

%%Plotting switch. 
switch flags.plot
    case 0
    %%Full output
    l=4; %grid height
    w=5; %grid width
    kk=1:20;
    case 1
    %%9 box output
    l=3; %grid height
    w=3; %grid width
    kk=[18,2,3,1,7,13,17,19,9];
    case 1.1
    %%9 box output
    l=3; %grid height
    w=3; %grid width
    kk=[18,2,3,1,7,13,21,19,9];
    case 2
    %%6 box output
    l=3; %grid height
    w=2; %grid width
    kk=[18,2,7,6,9,19];
    case 2.1
    %%6 box output
    l=2; %grid height
    w=3; %grid width
    kk=[18,2,7,19,6,9];
    case 2.2
    %%6 box output
    l=2; %grid height
    w=3; %grid width
    kk=[18,2,7,8,9,17];
    case 3
    %%4 box output
    l=2; %grid height
    w=2; %grid width
    kk=[6,2,13,7]; 
    case 3.1
    %%4 box output
    l=2; %grid height
    w=2; %grid width
    kk=[18,2,9,7];
    case 3.2
    %%4 box output
    l=2; %grid height
    w=2; %grid width
    kk=[19,2,9,7];
    case 3.3
    %%4 box output
    l=2; %grid height
    w=2; %grid width
    kk=[7,2,9,17];
    case 4
    %%4 box output
    l=1; %grid height
    w=3; %grid width
    kk=[18,2,9];    
    case 5
    %%2 box output
    l=1; %grid height
    w=2; %grid width
    kk=[18,2];
    case 5.1
    %%2 box output
    l=1; %grid height
    w=2; %grid width
    kk=[9,2];
    case 5.2
    %%2 box output
    l=2; %grid height
    w=1; %grid width
    kk=[7,9];
    case 5.3
    %%2 box output
    l=2; %grid height
    w=1; %grid width
    kk = 2;
    case 5.4
    %%2 box output
    l=3; %grid height
    w=1; %grid width
    kk = [18, 2] ;
    
    case 5.5
    %%2 box output
    l=3; %grid height
    w=1; %grid width
    kk = [7, 2] ;
    
    

end
    
    
%make some labels   
ylab={'DIC  [mmol/kg]'; '\delta^{13}C'; '[Ca^{2+}]   [mmol/kg]';
'F^b_{org} [mol/yr]';'F^b_{carb} [mol/yr]'; 
'f_{org}'; 'pCO_2  [ppmv]';'pH  ';  
'\Omega';  '\epsilon_b';  'HCO_3 [mol/kg]'; 
'Avg \delta^{13}C^w ';  'ALK [meq/kg]' ; 'Avg \delta^{13}C^b' 
'\delta^{13}C^b_{org}'; '[PO_4^{3-}]';'\delta^{44}Ca'...
;'F_{extra}  [mol/yr]'; 'C:P';'F^w_{sil} [mol/yr]';...
'[O_2] %'};





%lines denoting perturbation max (changes according to perturbation chosen)    
% switch flags.pertb
%    case 0
%         tforcemax=0;
%    case 1
%        tforcemax=Output.t(Output.Fborg==max(Output.Fborg));
%    case 2
%        tforcemax=Output.t(Output.Fwsil==max(Output.Fborg));
%    case {3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8}
%        tforcemax=Output.t(Output.Fextra==max(Output.Fextra));
% end  



switch flags.plot
    case {0,1,1.1,2,2.1,2.2,3,3.1,3.2,3.3,4,5,5.1,5.2,5.3,5.4,5.5}
        %plotting loop, makes subplots
        nn = length(kk); 
        for ii = 1:(nn)
                    hs(ii) = subplot(l,w,ii);
                     line(x,y(:,kk(ii)),'LineStyle','-');
                     %,'LineWidth',1.5,'Color','k')
                    xlabel('Time')
                     ylabel(ylab(kk(ii)))
                     set(gca,'XMinorTick','on')
                     grid on
                     box on
                     figure(hf)
                    
                
                    
%                      annotation('textbox',get(gca,'Position')-[0.05 -.06 0 0],'String',['(',char(ii+64),')']...
%                         ,'LineStyle','none', 'FitBoxToText','on');
%                      if flags.pertb~=3.9 && flags.plot~=4 && tcrit>0.2
%                       line([params.t_perturbation_start_1 params.t_perturbation_start_1],ylim,'Color','r','LineStyle','--')
%                       line([tforcemax tforcemax],ylim,'Color','g','LineStyle','-.')
%                       line([params.t_perturbation_end_1 params.t_perturbation_end_1],ylim,'Color','r','LineStyle','--')
%                      end
%                     axis tight 

        end



    case 6 

        a1 = gca;
        a2 = axes('Position',get(a1,'Position'),...
            'YAxisLocation','right','XAxisLocation','top','XTickLabel',...
            '','Color','none');
                
        yl1 = line(x,y(:,2),'Parent',a1);
        yl2 = line(x,y(:,15),'Parent',a2,'LineStyle','--','Color','k');
    
        ylabel(a1,'\delta^{13}C_{carb}')
        xlabel(a1,'time (yrs)')
        grid on 

        hl = legend([yl1 yl2],'\delta^{13}C_{carb}','\delta^{13}C_{org}');
        ylabel(a2,'\delta^{13}C_{org}')



    case 7
      hold on
      scatter(Output.DIC, Output.ALK,[],Output.t)   

end





%%save figure to directory as eps. 
switch flags.plot_save
    case 1    
    
        %%%Uncomment one or the other to use:
        
        %%print using screen size:   
        set(gcf,'PaperPositionMode','auto')
     
        %or resize figure to standard two column width (190mm) and a reasonable
        %height 
        
%         set(gcf, 'PaperUnits', 'centimeters')
%         papersize = get(gcf, 'PaperSize');
%         width = 19.0;         % Initialize a variable for width.
%         height = 15;          % Initialize a variable for height.
%         left = (papersize(1)- width)/2;
%         bottom = (papersize(2)- height)/2;
%         myfiguresize = [left, bottom, width, height];
%         set(gcf, 'PaperPosition', myfiguresize);
    
    %print figure 
    print('-f1','-depsc','-loose', '-painters', [flags.file_path flags.file_name]);
end




end

function [Output] = model(params,flags)
% this is the model function. Is called by the control function.  
% subfunction odefun holds the diff equations

   

%Oceanic parameters
tempi = 20; %[C] Celsius
depth = 0;  %[m] meters
sal = 1.035;   %[kg/L]
Voc = 1.32e21; %[L]

%initial mass of carbon

switch flags.pCO2
    
    case 0 
        DICi = 0.3e-3; %[mol/kg]
        ALKi = 0.9e-3; %[eq/kg]   
    case 1 
        DICi = 2.0e-3; %[mol/kg]
        ALKi = 2.2e-3; %[eq/kg]   
        
    case 2 
        DICi = 20.0e-3; %[mol/kg]
        ALKi = 19.0e-3; %[eq/kg] 
      
    case 3    
        DICi = 5.0e-3; %[mol/kg]
        ALKi = 5.22e-3; %[eq/kg] 
        
    case 4 
        DICi = 1.9e-3; %[mol/kg]
        ALKi = 1.89e-3; %[eq/kg] 
        
end   
   
       
MCi = DICi*sal*Voc; %mass of C in the ocean in [mol] 

%initial mass of calcium
Cai = 17e-3; %[mol/kg]  modern: 10.3e-3
% Cai = 14e-3; %[mol/kg]  modern: 10.3e-3
MCai=Cai*sal*Voc; %mass of Ca in the oceans in [mol] 


%initial alkalinity
kalk = 2*Cai - ALKi;%[eq/kg]
%kalk the net difference between the cations and anions 
%excluding Ca


%initial carbonate system parameters
[pco2i,pH,co2i,hco3,co3]=CO3eq(tempi,sal,depth,ALKi,DICi); %other params
omegai=(MCai/Voc/sal*co3)/(10^-6.37); %calcite saturation state. 

PO4i=0.25e-6 .* params.PO4_multi; %[mol/kg]



%photosynthetic fractionation
switch flags.del_diff
    case 0
    epsCi = 30;    %[permil]
    case 1
    epsCi = 30;
    case 2
    epsCi = -((159.5*PO4i+38.39)/(co2i*1e6) - 33);
    case 3
    epsCi = -((159.5*PO4i+38.39)/(co2i*1e6) - 33);
end


% initial Carbon input fluxes
fworg = 0.25;  
delvolc = -5.5;  %[permil]


Fvolc = 16e12; %[mol/y]
delwcarb = delvolc + fworg * epsCi; %[permil]
delworg = delwcarb-epsCi; %[permil]
Fwcarbi = 24e12; %[mol/y]
Fworgi = fworg*Fwcarbi/(1-fworg);  %[mol/y]


forgi = fworg;   
Fbcarbi = (1-forgi)*(Fvolc+Fworgi+Fwcarbi); %[mol/y]
Fborgi = (forgi)*(Fvolc+Fworgi+Fwcarbi);  %[mol/y]
Fborg = Fborgi;

%initial Ca input fluxes
Fhyd = 0; %[mol/y]
Fwsili = Fvolc + Fworgi - Fhyd - Fborgi;   %[mol/y]


%O2 fluxes
params.f_red = fworg;
MO2 = 1.8e20*0.2095;


%phosphate:
%initial mass of Phosphate
CPi=106;

Mpi=PO4i*sal*Voc; %mass of P in the ocean [mol]
Fwpi=(Fwcarbi+Fvolc+Fworgi-(Fhyd+Fwsili+Fwcarbi))/CPi; %[mol/yr]
Fbpi=Fwpi; %[mol/y]


epsCa=1.4;  %[permil]

%isotopic values of initial input fluxes
%carbon

%calcium
delhyd=-0.25;  %[permil]
delwriv=-0.6;  %[permil]


%isotopic value of initial output fluxes
delw=(Fvolc*delvolc+Fworgi*delworg+Fwcarbi*delwcarb)/(Fvolc+Fworgi+Fwcarbi);
delbcarbi = delw + forgi*epsCi;
delborgi = delbcarbi - epsCi;
delb = (Fbcarbi*delbcarbi+Fborgi*delborgi)/(Fborgi+Fbcarbi);
fborg_test = Fborgi/(Fbcarbi+Fborgi);
delbcarbCai = (delhyd*Fhyd + delwriv*(Fwsili + Fwcarbi))/...
    (Fhyd+Fwsili+Fwcarbi) + epsCa;

% disp([
% Mpi/Fbpi;
% MCi/Fborgi;
% MCi/Fbcarbi;
% MCai/Fbcarbi])


% Berner weathering params
    G=3.3;
    Z=0.09;
    Zc=0.07;

    


%%
%the perturbation


n=5000;     %number of points in the interpolation
tinterp=[params.t_run_start linspace(params.t_perturbation_start_1,...
    params.t_perturbation_end_1,n) params.t_run_end];


switch flags.pertb
case 0
    F_extra_interp= zeros(1,n+2);
case 1
    F_extra_interp=[Fborgi Fborgi*(1 + params.f_extra*normpdf(linspace(-4,4,n),0,1)/...
    max(normpdf(linspace(-4,4,n),0,1))) Fborgi];    
case 2
     F_extra_interp=[1 (1 + params.f_extra*normpdf(linspace(-4,4,n),0,1)/...
     max(normpdf(linspace(-4,4,n),0,1))) 1];    

case {3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.1,5,5.1}

    %first perturbation
    t_interp_1_pertb = ...
    linspace(params.t_perturbation_start_1,params.t_perturbation_end_1,n);

    t_interp_1 = [params.t_run_start t_interp_1_pertb(1) ...
    t_interp_1_pertb(2:(end -1)) t_interp_1_pertb(end)...
    params.t_run_end];

    
    sigma_1 = (params.t_perturbation_end_1 - params.t_perturbation_start_1)/10;
    %ten standard deviations out in each direction

    mu_1 = (params.t_perturbation_end_1 + params.t_perturbation_start_1)/2;

    F_extra_interp_1 = [0 0 params.M_extra_1*...
    normpdf(t_interp_1_pertb(2:(end -1)), mu_1, sigma_1) 0 0]; 

            
end


switch flags.pertb

case {3.9,4,4.1,5,5.1}


%second perturbation
t_interp_2_pertb = ...
    linspace(params.t_perturbation_start_2,params.t_perturbation_end_2,n);

t_interp_2 = [params.t_run_start t_interp_2_pertb(1) ...
    t_interp_2_pertb(2:(end -1)) t_interp_2_pertb(end)...
    params.t_run_end];

switch flags.run

case {'multi_pert','custom'}
sigma_2 = (params.t_perturbation_end_2 - params.t_perturbation_start_2)/4; 

mu_2 = (50*params.t_perturbation_start_2 + params.t_perturbation_end_2)/51 ;

a = 20;

F_extra_interp_2 = [0 0 params.M_extra_2*skew_norm(...
    t_interp_2_pertb(2:(end -1)), mu_2, sigma_2, a ) 0 0]; 



otherwise 

sigma_2 = (params.t_perturbation_end_2 - params.t_perturbation_start_2)/10; 
%ten standard deviations out in each direction

mu_2 = (params.t_perturbation_start_2 + params.t_perturbation_end_2)/2;
%             
F_extra_interp_2 = ...
[0 0 params.M_extra_2*normpdf(t_interp_2_pertb(2:(end -1)), mu_2, sigma_2) 0 0]; 

end

            
end

switch flags.pertb
                        
    case {4,4.1,5,5.1}    
    
        
    %third perturbation
    
    t_interp_3_pertb = ...
    linspace(params.t_perturbation_start_3,params.t_perturbation_end_3,n);

    t_interp_3 = [params.t_run_start t_interp_3_pertb(1) ...
    t_interp_3_pertb(2:(end -1)) t_interp_3_pertb(end)...
    params.t_run_end];

    
    switch flags.run
        
    case {'multi_pert','custom'}
        
        sigma_3 = (params.t_perturbation_end_3 - params.t_perturbation_start_3)/7; 

        mu_3 = (20*params.t_perturbation_start_3 + params.t_perturbation_end_3)/21;

        a = 4;

        F_extra_interp_3 = [0 0 params.M_extra_3*skew_norm( ...
            t_interp_3_pertb(2:(end -1)), mu_3, sigma_3, a ) 0 0]; 



    otherwise 

    sigma_3 = (params.t_perturbation_end_3 - params.t_perturbation_start_3)/10; 
    %ten standard deviations out in each direction

    mu_3 = (params.t_perturbation_start_3 + params.t_perturbation_end_3)/2;
    %             
    F_extra_interp_3 = ...
    [0 0 params.M_extra_3*normpdf(t_interp_3_pertb(2:(end -1)), mu_3, sigma_3) 0 0]; 

    end
    

    
    
    
    
    %fourth perturbation
    
    t_interp_4_pertb = ...
    linspace(params.t_perturbation_start_4,params.t_perturbation_end_4,n);

    t_interp_4 = [params.t_run_start t_interp_4_pertb(1) ...
    t_interp_4_pertb(2:(end -1)) t_interp_4_pertb(end)...
    params.t_run_end];

    
    switch flags.run
        
    case {'multi_pert','custom'}
        
        sigma_4 = (params.t_perturbation_end_4 - params.t_perturbation_start_4)/10; 

        mu_4 = (20*params.t_perturbation_start_4 + params.t_perturbation_end_4)/21;

        a = 4;

        F_extra_interp_4 = [0 0 params.M_extra_4*skew_norm( ...
            t_interp_4_pertb(2:(end -1)), mu_4, sigma_4, a ) 0 0]; 



    otherwise 

    sigma_4 = (params.t_perturbation_end_4 - params.t_perturbation_start_4)/10; 
    %ten standard deviations out in each direction

    mu_4 = (params.t_perturbation_start_4 + params.t_perturbation_end_4)/2;
    %             
    F_extra_interp_4 = ...
    [0 0 params.M_extra_4*normpdf(t_interp_4_pertb(2:(end -1)), mu_4, sigma_4) 0 0]; 

    end

            
end


switch flags.pertb
    case 4.1
        %silicate weathering shift
            
            
            t_shift = linspace(params.t_shift_start,params.t_run_end,n);
            
            t_interp_shift = [params.t_run_start  t_shift ];

            Fwsil_extra_interp = [0 0 params.shift_mag*ones(1,length(t_shift)-1)];
end

switch flags.pertb
            
     case {5,5.1}
            
         %fifth perturbation
            t_interp_5_pertb = ...
            linspace(params.t_perturbation_start_5,params.t_perturbation_end_5,n); 
            
            t_interp_5 = [params.t_run_start t_interp_5_pertb(1) ...
            t_interp_5_pertb(2:(end -1)) t_interp_5_pertb(end)...
            params.t_run_end];
 
            
            switch flags.run
            
            case {'custom','multi_pert'}
                
            
            sigma_5 = (params.t_perturbation_end_5 - params.t_perturbation_start_5)/10; 
            %ten standard deviations out in each direction
            
            mu_5 = (20*params.t_perturbation_start_5 + params.t_perturbation_end_5)/21;
                             
            a = 4;
            
            F_extra_interp_5 = [0 0 params.M_extra_5*skew_norm(...
                t_interp_5_pertb(2:(end -1)), mu_5, sigma_5, a ) 0 0]; 
            
            otherwise 
            
                       
            sigma_5 = (params.t_perturbation_end_5 - params.t_perturbation_start_5)/8; 
            %ten standard deviations out in each direction
            
            mu_5 = (params.t_perturbation_end_5 + params.t_perturbation_start_5)/2;
            
            F_extra_interp_5 = [0 0 ...
                params.M_extra_5*normpdf(t_interp_5_pertb(2:(end -1)),...
                mu_5, sigma_5) 0 0]; 
                        
            end
         
         
            
        
            %silicate weathering shift
            
            
            t_shift = linspace(params.t_shift_start,params.t_run_end,n);
            
            t_interp_shift = [params.t_run_start  t_shift ];

            Fwsil_extra_interp = [0 0 params.shift_mag*ones(1,length(t_shift)-1)];
            
            
            %CP power shift
            switch flags.pertb 
                case 5.1
                t_interp_CP_bP = [params.t_run_start,...
                 params.t_p_CP_bP - 1e3, ...
                 params.t_p_CP_bP + 1e3, params.t_run_end];
        
                F_extra_interp_CP_bP = ...
              [ params.p_CP_bP params.p_CP_bP ... 
                params.p_CP_bP_2 params.p_CP_bP_2 ];
            
                        
                        
            
            end
    

end                   

%%

%the solver
ic=[MCi delbcarbi MCai delbcarbCai Mpi MO2]; %initial conditions
% options=odeset('RelTol',1e-7,'OutputFcn',@stop_func);


switch flags.pertb

    case 3.9
    
    [t1 y1] = ode15s(@odefun,[params.t_run_start params.t_perturbation_start_1],ic); %the solver up to the pertubation
      
    options_2 = odeset('InitialStep', mean(diff(t_interp_1_pertb)));
    
    [t2 y2] = ode15s(@odefun, [params.t_perturbation_start_1 params.t_perturbation_start_2], y1(end,:), options_2); 
    
    options_3 = odeset('InitialStep', mean(diff(t_interp_2_pertb)) );
    
    [t3 y3] = ode15s(@odefun, [params.t_perturbation_start_2 params.t_run_end], y2(end,:), options_3);  
    
    
    %concatenate the solutions
    t = [t1; t2; t3];
    y = [y1; y2; y3 ];
    
    
    case {4,4.1}
    
    [t1 y1] = ode15s(@odefun,[params.t_run_start params.t_perturbation_start_1],ic); %the solver up to the pertubation
      
    options_2 = odeset('InitialStep', mean(diff(t_interp_1_pertb)) );
    
    [t2 y2] = ode15s(@odefun, [params.t_perturbation_start_1 params.t_perturbation_start_2], y1(end,:), options_2); 
    
    options_3 = odeset('InitialStep', mean(diff(t_interp_2_pertb)) );
    
    [t3 y3] = ode15s(@odefun, [params.t_perturbation_start_2 params.t_perturbation_start_3], y2(end,:), options_3); 
    
    options_4 = odeset('InitialStep', mean(diff(t_interp_3_pertb)) );
    
    [t4 y4] = ode15s(@odefun, [params.t_perturbation_start_3 params.t_perturbation_start_4], y3(end,:), options_4); 
    
    options_5 = odeset('InitialStep', mean(diff(t_interp_4_pertb)) );
    
    [t5 y5] = ode15s(@odefun, [params.t_perturbation_start_4 params.t_run_end], y4(end,:), options_5);  
    
    
    %concatenate the solutions
    t = [t1; t2; t3; t4; t5];
    y = [y1; y2; y3; y4; y5];    
        
    
    case 5
    
    [t1 y1] = ode15s(@odefun,[params.t_run_start params.t_perturbation_start_1],ic); %the solver up to the pertubation
    
    
    options_2 = odeset('InitialStep', mean(diff(t_interp_1_pertb)) );
   
    [t2 y2] = ode15s(@odefun, [params.t_perturbation_start_1 params.t_perturbation_start_2], y1(end,:), options_2); 
    
    
    options_3 = odeset('InitialStep', mean(diff(t_interp_2_pertb)) );
    
    [t3 y3] = ode15s(@odefun, [params.t_perturbation_start_2 params.t_perturbation_start_3], y2(end,:), options_3); 
    
    
    options_4 = odeset('InitialStep', mean(diff(t_interp_3_pertb)) );
    
    [t4 y4] = ode15s(@odefun, [params.t_perturbation_start_3 params.t_perturbation_start_4], y3(end,:), options_4); 
    
    
    options_5 = odeset('InitialStep', mean(diff(t_interp_4_pertb)) );
    
    [t5 y5] = ode15s(@odefun, [params.t_perturbation_start_4 params.t_perturbation_start_5], y4(end,:), options_5);  
    
    
    options_6 = odeset('InitialStep', mean(diff(t_interp_5_pertb)) );
    
    [t6 y6] = ode15s(@odefun, [params.t_perturbation_start_5 params.t_run_end], y5(end,:), options_6);  
    
    
    %concatenate the solutions
    t = [t1; t2; t3; t4; t5; t6];
    y = [y1; y2; y3; y4; y5; y6];
    
    
    
    otherwise
    min_duration = min(params.t_perturbation_end_1 - params.t_perturbation_start_1);    
    options1 = odeset('Refine',4);
    options2 = odeset('InitialStep',min_duration/10,'Refine',4);
        
    [t1 y1]=ode15s(@odefun,[params.t_run_start params.t_perturbation_start_1],ic,options1); %the solver up to the pertubation
    [t2 y2]=ode15s(@odefun,[params.t_perturbation_start_1 params.t_run_end],...
        y1(end,:), options2 ); %the solver after the perturbation   
    
    %concatenate the two solutions
    t=[t1; t2];
    y=[y1; y2];
       
    
end




%%
%the solver calls this function 
    function dy=odefun(t,y)
    

    %carbonate system
    ALK=(2*(y(3)/(sal*Voc))-kalk); %[eq/kg] carbonate alkalinity
    DIC=y(1)/(sal*Voc); %[mol/kg] dissolved inorganic carbon
    [pco2,pH,co2,hco3,co3]=CO3eq(tempi ,sal,depth,ALK,DIC) ;
    %omega=(y(3)/Voc/sal*1e3*co3)/(10^-6.37) %calcite saturation state. 
    omega=(y(3)/Voc/sal*co3)/(10^-6.37);
    RCO2=pco2/pco2i;
    
    
    %fractionation dependancies
    switch flags.del_diff
        case 0
            epsC = 30;    %[permil]
        case 1
            epsC = 30;
        case 2
            epsC = -((159.5*PO4i+38.39)/(co2*1e6) - 33) ;
        case 3
            epsC = -((159.5*y(5)/(sal*Voc)+38.39)/(co2*1e6) - 33) ;
    end
    
    

     %weathering flux dependancies
    switch flags.weathering
        case 0
            Fwcarb=Fwcarbi;
            Fworg=Fworgi;
            Fwsil=Fwsili;
            Fwp=Fwpi;
       case 1
           Fwcarb=Fwcarbi*(RCO2^params.p_w);
           Fworg=Fworgi*(RCO2^params.p_w);
           Fwsil=Fwsili*(RCO2^params.p_w);
           Fwp=Fwpi*(RCO2^params.p_w);
            
        case 2
            Fwcarb=Fwcarbi*(RCO2)^(G*Z)*(1+G*Zc*log(RCO2));
            Fworg=Fworgi*(RCO2)^(G*Z)*(1+G*Z*log(RCO2))^0.65;
            Fwsil=Fwsili*(RCO2)^(G*Z)*(1+G*Z*log(RCO2))^0.65;
            Fwp=Fwpi*(RCO2)^(G*Z)*(1+G*Z*log(RCO2))^0.65;
    end
   
    
    
            
    
    
   %phophate burial 
   switch flags.PO4
       case 0
            Fbp = Fbpi * (y(5) / Mpi).^params.p_bp_Mp;
   end
   
   %the perturbation (organic carbon burial and extra carbon)
    
   Fwsil_extra = 0;
   
   switch flags.pertb        
        case 0
            F_extra = 0;
            del_extra = params.del_extra_1;
            Fbcarb=Fbcarbi;
            Fborg=Fborgi;
            
        case 1
            F_extra = 0;
            del_extra = params.del_extra_1;
            Fborg = interp1(tinterp,F_extra_interp,t,'pchip');
            
        case 2
            F_extra = 0;
            del_extra = params.del_extra_1;
            CP = CPi * (Fbp/Fbpi).^params.p_CP_bP;
            Fborg = Fbp*CP; %[mol/y]
            Fwsil = Fwsili*interp1(tinterp,F_extra_interp,t,'pchip');
            Fwcarb=Fwcarbi*interp1(tinterp,F_extra_interp,t,'pchip');
            Fworg=Fworgi*interp1(tinterp,F_extra_interp,t,'pchip');
            Fwp=Fwpi*interp1(tinterp,F_extra_interp,t,'pchip');
                             
        case 3.1
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            del_extra = params.del_extra_1;
            Fborg = Fborgi * RCO2; %[mol/y]
        case 3.2
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            del_extra = params.del_extra_1;
            CP = CPi;
            Fborg = CP*Fbp; %[mol/y]
        case 3.3
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            del_extra = params.del_extra_1;
%             CP = CPi * (RCO2)^p_CP_RCO2;
%             CPi * log (Fborg/Fborgi)
            CP = CPi + 0.05 * log (Fborg/Fborgi); 
            Fborg = CP*Fbp; %[mol/y]
       case 3.4 
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            del_extra = params.del_extra_1;
            CP = CPi;
            Fborg = Fborgi; %[mol/y]
       
       case 3.6
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            del_extra = params.del_extra_1;
            CP = CPi * (Fbp/Fbpi).^params.p_CP_bP;
            Fborg = Fbp*CP; %[mol/y]
            
            
       case 3.7
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            del_extra = params.del_extra_1;
            CP = CPi * exp(Fbp / Fbpi -1 );
            Fborg = Fbp*CP; %[mol/y]
       case 3.8
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            del_extra = params.del_extra_1;
            CP = CPi;
            Fborg = Fwp.*CP; %[mol/y]
       case 3.9
            F_extra_1 = interp1(t_interp_1, F_extra_interp_1,t,'pchip');
            F_extra_2 = interp1(t_interp_2, F_extra_interp_2,t,'pchip');
            F_extra = F_extra_1 + F_extra_2;
            
            del_extra = (F_extra_1*params.del_extra_1 + F_extra_2*params.del_extra_2)/(F_extra);
            if isnan(del_extra) == 1, del_extra = 0; end
            
            CP =  CPi * (Fbp/Fbpi).^params.p_CP_bP;
            Fborg = Fbp*CP; %[mol/y]

            
        case 4
            F_extra_1 = interp1(t_interp_1, F_extra_interp_1,t,'pchip');
            F_extra_2 = interp1(t_interp_2, F_extra_interp_2,t,'pchip');
            F_extra_3 = interp1(t_interp_3, F_extra_interp_3,t,'pchip');
            F_extra_4 = interp1(t_interp_4, F_extra_interp_4,t,'pchip');
            
            F_extra = F_extra_1 + F_extra_2 + F_extra_3 + F_extra_4;
            
            del_extra = ...
            (F_extra_1*params.del_extra_1 + F_extra_2*params.del_extra_2 +...
            F_extra_3*params.del_extra_3 + F_extra_4*params.del_extra_4)...
            /(F_extra);
            
            if isnan(del_extra) == 1, del_extra = 0; end
            
            CP =  CPi * (Fbp/Fbpi).^params.p_CP_bP;
            Fborg = Fbp*CP; %[mol/y]
            
       case 4.1
           
            F_extra_1 = interp1(t_interp_1, F_extra_interp_1,t,'pchip');
            F_extra_2 = interp1(t_interp_2, F_extra_interp_2,t,'pchip');
            F_extra_3 = interp1(t_interp_3, F_extra_interp_3,t,'pchip');
            F_extra_4 = interp1(t_interp_4, F_extra_interp_4,t,'pchip');
            
            
            F_extra = F_extra_1 + F_extra_2 +...
                F_extra_3 + F_extra_4;
            
            del_extra = ...
            (F_extra_1*params.del_extra_1 + F_extra_2*params.del_extra_2 +...
            F_extra_3*params.del_extra_3 + F_extra_4*params.del_extra_4)...
            /(F_extra);
            
            if isnan(del_extra) == 1, del_extra = 0; end
            
            CP =  CPi * (Fbp/Fbpi).^params.p_CP_bP;
            Fborg = Fbp*CP; %[mol/y]
       
            Fwsil_extra = interp1(t_interp_shift, Fwsil_extra_interp,...
                t,'pchip');
            
         case 5
            F_extra_1 = interp1(t_interp_1, F_extra_interp_1,t,'pchip');
            F_extra_2 = interp1(t_interp_2, F_extra_interp_2,t,'pchip');
            F_extra_3 = interp1(t_interp_3, F_extra_interp_3,t,'pchip');
            F_extra_4 = interp1(t_interp_4, F_extra_interp_4,t,'pchip');
            F_extra_5 = interp1(t_interp_5, F_extra_interp_5,t,'pchip');
            
            F_extra = F_extra_1 + F_extra_2 +...
                F_extra_3 + F_extra_4 + F_extra_5;
            
            del_extra = ...
            (F_extra_1*params.del_extra_1 + F_extra_2*params.del_extra_2 +...
            F_extra_3*params.del_extra_3 + F_extra_4*params.del_extra_4...
            + F_extra_5*params.del_extra_5)...
            /(F_extra);
            
            if isnan(del_extra) == 1, del_extra = 0; end
            
            CP =  CPi * (Fbp/Fbpi).^params.p_CP_bP;
            Fborg = Fbp*CP; %[mol/y]
       
            Fwsil_extra = interp1(t_interp_shift, Fwsil_extra_interp,...
                t,'pchip');
            
         case 5.1
            F_extra_1 = interp1(t_interp_1, F_extra_interp_1,t,'pchip');
            F_extra_2 = interp1(t_interp_2, F_extra_interp_2,t,'pchip');
            F_extra_3 = interp1(t_interp_3, F_extra_interp_3,t,'pchip');
            F_extra_4 = interp1(t_interp_4, F_extra_interp_4,t,'pchip');
            F_extra_5 = interp1(t_interp_5, F_extra_interp_5,t,'pchip');
            
            F_extra = F_extra_1 + F_extra_2 +...
                F_extra_3 + F_extra_4 + F_extra_5;
            
            del_extra = ...
            (F_extra_1*params.del_extra_1 + F_extra_2*params.del_extra_2 +...
            F_extra_3*params.del_extra_3 + F_extra_4*params.del_extra_4...
            + F_extra_5*params.del_extra_5)...
            /(F_extra);
            
            if isnan(del_extra) == 1, del_extra = 0; end
            
            p_CP_bP = interp1(t_interp_CP_bP, F_extra_interp_CP_bP,t);
            
            CP =  CPi * (Fbp/Fbpi).^ p_CP_bP;
            Fborg = Fbp*CP; %[mol/y]
       
            Fwsil_extra = interp1(t_interp_shift, Fwsil_extra_interp,...
                t,'pchip');   
            
    end
    
    
    %carbonate buiral
    switch flags.Fbcarb
        case 0
          Fbcarb = Fbcarbi;
        case 1
          Fbcarb = Fhyd + Fwsil  + Fwcarb ;
        case 2
          Fbcarb = Fwcarb + Fworg + Fvolc + F_extra - Fborg;
        case 3
          Fbcarb = Fbcarbi * y(3)/MCai; %[mol/y]  
        case 4
          Fbcarb = Fbcarbi*(omega/omegai).^params.p_carb;
       case 5
          Fbcarb = Fbcarbi*( (omega-1)/(omegai-1) ).^params.p_carb;
          
    end
    
    
    switch flags.pertb 
%         case 2
%             Fbcarb = interp1(tinterp,F_extra_interp,t,'pchip');
%             Fborg = Fwcarb + Fworg + Fvolc + F_extra - Fbcarb;
        case 3.0
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            Fborg = Fwcarb + Fworg + Fvolc + F_extra - Fbcarb; %[mol/y] 
        case 3.5 
           F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            CP = CPi;
            Fborg = forgi/(1-forgi)*Fbcarb; %[mol/y]
            del_extra = params.del_extra_1;
    end 
    
    

    dy(1) = Fwcarb + Fworg + Fvolc + F_extra - Fborg - Fbcarb; %d(MC)/dt
    dy(2) = ( Fvolc*(delvolc-y(2)) + Fworg*(delworg-y(2)) +...
        Fwcarb*(delwcarb-y(2)) + F_extra*(del_extra - y(2))...
        - (Fborg)*(-epsC)) / y(1); %d(delC)/dt
    
    
    %Calcium boxes
    dy(3) = Fhyd + Fwsil + Fwsil_extra + Fwcarb - Fbcarb; %d(MCa)/dt
    dy(4) = ( Fhyd*(delhyd-y(4)) + (Fwsil + Fwsil_extra + Fwcarb)*(delwriv-y(4))...
        - Fbcarb*(-epsCa) ) / y(3); %d(delCa)/dt
    
    %PO4 box
    dy(5) = Fwp - Fbp;
    
    %O2 box
    dy(6) = Fborg - Fworg - params.f_red .* Fvolc;
   
  

    dy=dy(:);
  
    
    end


    
%%
%output

%reconstruct vars
MC = y(:,1);
MCa = y(:,3);
DIC = y(:,1)/(sal*Voc); %[mmol/kg] dissolved inorganic carbon
Ca = y(:,3)/(sal*Voc);
delbcarb = y(:,2);
delCa = y(:,4);
PO4 = y(:,5)/(sal*Voc);
Mp = y(:,5);
MO2 = y(:,6);



ALK=(2*Ca-kalk); %[meq/kg] carbonate alk
pco2=zeros(length(ALK),1);pH=zeros(length(ALK),1);omega=zeros(length(ALK),1);
co2=zeros(length(ALK),1); hco3=zeros(length(ALK),1); co3=zeros(length(ALK),1);
for i=1:length(ALK),
    [pco2(i),pH(i),co2(i),hco3(i),co3(i)]=CO3eq(tempi,sal,depth,ALK(i),DIC(i)); 
    omega(i)=(Ca(i).*co3(i))/(10^-6.37); %calcite saturation state. 
    
end

    switch flags.del_diff
        case 0
        epsC = 30;    %[permil]
        delborg = delbcarb - epsC;
        case 1
        epsC = 30;
        delborg = (1-params.fexternal)*(delbcarb - epsC) + params.fexternal*params.delexternal;
        case 2
        epsC = -((159.5.*PO4i+38.39)./(co2*1e6) - 33) ;
        delborg = delbcarb - epsC;
        case 3
        epsC = -((159.5.*PO4+38.39)./(co2*1e6) - 33) ;
        delborg = delbcarb - epsC;
    end

    
    
    RCO2=pco2/pco2i;
    
    %weathering feedbacks
    switch flags.weathering
            
        case 0
            Fwcarb=Fwcarbi*ones(length(t),1);
            Fworg=Fworgi*ones(length(t),1);
            Fwsil=Fwsili*ones(length(t),1);
            Fwp=Fwpi*ones(length(t),1);
        case 1
            Fwcarb=Fwcarbi.*(RCO2.^params.p_w);
            Fworg=Fworgi.*(RCO2.^params.p_w);
            Fwsil=Fwsili.*(RCO2.^params.p_w);
            Fwp=Fwpi.*(RCO2.^params.p_w);
        case 2
            Fwcarb=Fwcarbi.*(RCO2).^(G.*Z).*(1+G.*Zc.*log(RCO2));
            Fworg=Fworgi.*(RCO2).^(G.*Z).*(1+G.*Z.*log(RCO2)).^0.65;
            Fwsil=Fwsili.*(RCO2).^(G.*Z).*(1+G.*Z.*log(RCO2)).^0.65;
            Fwp=Fwpi.*(RCO2).^(G.*Z).*(1+G.*Z.*log(RCO2)).^0.65;
    end
       
    
    %phosphate burial
    switch flags.PO4
       case 0
            Fbp = Fbpi .* (Mp ./ Mpi).^params.p_bp_Mp;
    end
    
    
    
    %forcing
    switch flags.pertb        
        case 0
            F_extra = zeros(length(t),1);
            del_extra=zeros(length(t),1);
            Fborg=Fborgi*ones(length(t),1);
            CP = CPi * ones(length(t),1);
        case 1
            F_extra = zeros(length(t),1);
            del_extra=zeros(length(t),1);
            Fborg = interp1(tinterp,F_extra_interp,t,'pchip');
            CP = CPi * ones(length(t),1);
        case 2
            F_extra = zeros(length(t),1);
            del_extra=zeros(length(t),1);
            CP = CPi .* (Fbp./Fbpi).^params.p_CP_bP;
            Fborg = Fbp.*CP; %[mol/y]
            Fwsil = Fwsili*interp1(tinterp,F_extra_interp,t,'pchip');
            Fwcarb=Fwcarbi*interp1(tinterp,F_extra_interp,t,'pchip');
            Fworg=Fworgi*interp1(tinterp,F_extra_interp,t,'pchip');
            Fwp=Fwpi*interp1(tinterp,F_extra_interp,t,'pchip');
        case 3.0
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            CP = CPi * ones(length(t),1);
        case 3.1
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            Fborg = Fborgi * MC / MCi; %[mol/y]
            CP = CPi * ones(length(t),1);    
        case 3.2
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            CP = CPi * ones(length(t),1);
            Fborg = CP.*Fbp; %[mol/y]
            
        case 3.3
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
%             CP = CPi * (RCO2).^p_CP_RCO2;
%             CP = CPi + CPi .* log (Fborg./Fborgi) ;
%             Fborg = CP*Fbp; %[mol/y];
            Fborg = CP.*Fbp; %[mol/y]  
            
        case 3.4 
           F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
           CP = CPi * ones(length(t),1);
           Fborg = Fborgi * ones(length(t),1); %[mol/y]
        case 3.6
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            CP = CPi .* (Fbp./Fbpi).^params.p_CP_bP;
            Fborg = Fbp.*CP; %[mol/y]
        case 3.7
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            CP = CPi .* exp(Fbp ./ Fbpi -1 );
            Fborg = Fbp.*CP; %[mol/y]
            
        case 3.8
            F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            CP = CPi * ones(length(t),1);
            Fborg = Fwp.*CP; %[mol/y]
        
        case 3.9
            F_extra_1 = interp1(t_interp_1, F_extra_interp_1,t,'pchip');
            F_extra_2 = interp1(t_interp_2, F_extra_interp_2,t,'pchip');
            F_extra = F_extra_1 + F_extra_2;
            
            del_extra = (F_extra_1*params.del_extra_1 + F_extra_2*params.del_extra_2)./(F_extra);
            if isnan(del_extra) == 1, del_extra = 0; end
            
            CP =  CPi * (Fbp/Fbpi).^params.p_CP_bP;
            Fborg = Fbp.*CP; %[mol/y]
        case 4
            F_extra_1 = interp1(t_interp_1, F_extra_interp_1,t,'pchip');
            F_extra_2 = interp1(t_interp_2, F_extra_interp_2,t,'pchip');
            F_extra_3 = interp1(t_interp_3, F_extra_interp_3,t,'pchip');
            F_extra_4 = interp1(t_interp_4, F_extra_interp_4,t,'pchip');
            
            F_extra = F_extra_1 + F_extra_2 + F_extra_3 + F_extra_4;
            
            del_extra = ...
            (F_extra_1*params.del_extra_1 + F_extra_2*params.del_extra_2 +...
            F_extra_3*params.del_extra_3 + F_extra_4*params.del_extra_4)...
            ./(F_extra);
            
            if isnan(del_extra) == 1, del_extra = 0; end
            
            CP =  CPi * (Fbp/Fbpi).^params.p_CP_bP;
            Fborg = Fbp.*CP; %[mol/y]
        
        case 4.1
            
            F_extra_1 = interp1(t_interp_1, F_extra_interp_1,t,'pchip');
            F_extra_2 = interp1(t_interp_2, F_extra_interp_2,t,'pchip');
            F_extra_3 = interp1(t_interp_3, F_extra_interp_3,t,'pchip');
            F_extra_4 = interp1(t_interp_4, F_extra_interp_4,t,'pchip');
            
            F_extra = F_extra_1 + F_extra_2 + F_extra_3 + F_extra_4;
            
            del_extra = ...
            (F_extra_1*params.del_extra_1 + F_extra_2*params.del_extra_2 +...
            F_extra_3*params.del_extra_3 + F_extra_4*params.del_extra_4)...
            ./(F_extra);
            
            del_extra(isnan(del_extra)) = 0;
            
            
            CP =  CPi * (Fbp/Fbpi).^params.p_CP_bP;
            Fborg = Fbp.*CP; %[mol/y]    
    
            Fwsil_extra = interp1(t_interp_shift, Fwsil_extra_interp,...
                t,'pchip');
            
        case 5
            F_extra_1 = interp1(t_interp_1, F_extra_interp_1,t,'pchip');
            F_extra_2 = interp1(t_interp_2, F_extra_interp_2,t,'pchip');
            F_extra_3 = interp1(t_interp_3, F_extra_interp_3,t,'pchip');
            F_extra_4 = interp1(t_interp_4, F_extra_interp_4,t,'pchip');
            F_extra_5 = interp1(t_interp_5, F_extra_interp_5,t,'pchip');
            
            F_extra = F_extra_1 + F_extra_2 + F_extra_3 + F_extra_4...
                + F_extra_5;
            
            del_extra = ...
            (F_extra_1*params.del_extra_1 + F_extra_2*params.del_extra_2 +...
            F_extra_3*params.del_extra_3 + F_extra_4*params.del_extra_4...
            + F_extra_5*params.del_extra_5)...
            ./(F_extra);
            
            del_extra(isnan(del_extra)) = 0;
            
            
            CP =  CPi * (Fbp/Fbpi).^params.p_CP_bP;
            Fborg = Fbp.*CP; %[mol/y]    
    
            Fwsil_extra = interp1(t_interp_shift, Fwsil_extra_interp,...
                t,'pchip');
    
        case 5.1
            F_extra_1 = interp1(t_interp_1, F_extra_interp_1,t,'pchip');
            F_extra_2 = interp1(t_interp_2, F_extra_interp_2,t,'pchip');
            F_extra_3 = interp1(t_interp_3, F_extra_interp_3,t,'pchip');
            F_extra_4 = interp1(t_interp_4, F_extra_interp_4,t,'pchip');
            F_extra_5 = interp1(t_interp_5, F_extra_interp_5,t,'pchip');
            
            F_extra = F_extra_1 + F_extra_2 +...
                F_extra_3 + F_extra_4 + F_extra_5;
            
            del_extra = ...
            (F_extra_1*params.del_extra_1 + F_extra_2*params.del_extra_2 +...
            F_extra_3*params.del_extra_3 + F_extra_4*params.del_extra_4...
            + F_extra_5*params.del_extra_5)...
            ./(F_extra);
            
            if isnan(del_extra) == 1, del_extra = 0; end
            
            p_CP_bP = interp1(t_interp_CP_bP, F_extra_interp_CP_bP,t);
            
            CP =  CPi .* (Fbp/Fbpi).^ p_CP_bP;
            Fborg = Fbp.*CP; %[mol/y]
       
            Fwsil_extra = interp1(t_interp_shift, Fwsil_extra_interp,...
                t,'pchip'); 
            
            
    end    
    
 

    
    %carbonate burial
    switch flags.Fbcarb
        case 0
          Fbcarb = Fbcarbi* ones(length(t),1);
        case 1
          Fbcarb = Fhyd + Fwsil  + Fwcarb ;
        case 2
          Fbcarb = Fwcarb + Fworg + Fvolc + F_extra - Fborg;
        case 3
          Fbcarb = Fbcarbi * MCa/MCai; %[mol/y]   
        case 4
          Fbcarb = Fbcarbi.*(omega./omegai).^params.p_carb;
        case 5
          Fbcarb = Fbcarbi*( (omega-1)/(omegai-1) ).^params.p_carb;
    end
    
    
    switch flags.pertb
%         case 2
%             Fbcarb = interp1(tinterp,F_extra_interp,t,'pchip');
%             Fborg = Fwcarb + Fworg + Fvolc + F_extra - Fbcarb;
        case 3.0
            Fborg = Fwcarb + Fworg + Fvolc + F_extra - Fbcarb; %[mol/y]
        case 3.5 
           F_extra = interp1(tinterp, F_extra_interp_1,t,'pchip');
            CP = CPi* ones(length(t),1);
            Fborg = forgi./(1-forgi).*Fbcarb; %[mol/y]
    end
    
CO2 = MO2/1.8e20;    

delwC=(Fvolc*delvolc+Fworg*delworg+Fwcarb*delwcarb+F_extra.*del_extra)./(Fvolc+Fworg+Fwcarb+F_extra);

delbC = (Fbcarb.*delbcarb+Fborg.*delborg)./(Fborg+Fbcarb);

epsb = delwC - delbC;

forg=Fborg./(Fborg+Fbcarb);

Output = struct('t', t, 'DIC',DIC*1e3, 'delbcarb', delbcarb, 'Ca', Ca*1e3,...
    'Fborg', Fborg, 'Fbcarb', Fbcarb, 'forg', forg,...
    'pco2', pco2, 'pH', pH,  'omega', omega, 'epsb',  epsb,...
    'hco3', hco3, 'co3', co3, 'ALK', ALK*1e3, ...
    'delbC', delbC, 'delborg', delborg,'PO4',PO4,...
    'delCa',delCa,'Fextra',F_extra,'CP',CP,'Fwsil',Fwsil,'delwC',delwC,...
    'CO2',CO2);


function [F] = skew_norm(x,mu,sigma,a)

t = (x - mu) / sigma;

F = 2 / sigma * normpdf(t,0,1) .* normcdf(a*t,0,1);


end






 
    end    

function [ pco2,pH,co2,hco3,co3 ] = CO3eq( temp,s,z,alk,dic )
%CO3eq calculates carbonate system parameters from alkalinity and DIC. 
%Modified from Chemical Oceanography and the Marine Carbon Cycle 2008 by Emerson
%and Hedges p129 and Richard Zeebe's csys code. 
%%
%inputs temp degC sal in ppt depth in m alk in eq/kg dic in mol/kg
%outputs co3 hco3 co2 in mol/kg


% Test parameters.
% Simply uncomment and run in cell mode (ctrl+enter) for outputs

% temp=20; %[degC] 
% s=1.035;%[kg/L]
% z=0; %[m]
% alk=2.2e-3 ;%[eq/kg]
% dic=2.0e-3 ;%[mol/kg]

%params
s=s*1e3-1e3;
t=temp+273.15; 
Pr=z/10;     
R=83.131;
tbor=0.000416*s/35.0;

%constants
U1=-60.2409+93.4517*(100/t)+23.3585*log(t/100);
U2=s*(0.023517-0.023656*(t/100)+0.0047036*(t/100)^2);
KH=exp(U1+U2);
KB=exp((-8966.9-2890.53*s^0.5-77.942*s+1.728*s^1.5-0.0996*s^2)/t...
    +148.0248+137.1942*s^0.5+1.62142*s-(24.4344+25.085*s^0.5+...
    0.2474*s)*log(t)+0.053105*s^0.5*t);
K1=10^(-(3633.86/t-61.2172+9.67770*log(t)-0.011555*s+0.0001152*s^2));
K2=10^(-(471.78/t+25.9290-3.16967*log(t)-0.01781*s+0.0001122*s^2));
dvB= -29.48+0.1622*temp-0.002608*(temp)^2;
dv1=-25.20+0.1271*temp;
dv2=-15.82-0.0219*temp;
dkB=-0.00284;
dk1=-0.00308+0.0000877*temp;
dk2=0.00113-0.0001475*temp;
KB=(exp(-(dvB/(R*t))*Pr+(0.5*dkB/(R*t))*Pr^2))*KB;
K1=(exp(-(dv1/(R*t))*Pr+(0.5*dk1/(R*t))*Pr^2))*K1;
K2=(exp(-(dv2/(R*t))*Pr+(0.5*dk2/(R*t))*Pr^2))*K2;
KW1=148.96502-13847.26/t-23.65218*log(t);
KW2=(118.67/t-5.977+1.0495*log(t))*s^0.5-0.01615*s;
KW=exp(KW1+KW2);

%solve for H ion
a1=1;
a2=(alk+KB+K1);
a3=(alk*KB-KB*tbor-KW+alk*K1+K1*KB+K1*K2-dic*K1);
a4=(-KW*KB+alk*KB*K1-KB*tbor*K1-KW*K1+alk*K1*K2+KB*K1*K2-dic*KB*K1-2*dic*K1*K2);
a5=(-KW*KB*K1+alk*KB*K1*K2-KW*K1*K2-KB*tbor*K1*K2-2*dic*KB*K1*K2);
a6=-KB*KW*K1*K2;
p=[a1 a2 a3 a4 a5 a6];
r=roots(p);
h=max(real(r));

%Calculate HCO3, CO3 ,CO2aq [mol/kg], and pCO2 [microatm]  
hco3=dic/(1+h/K1+K2/h);
co3=dic/(1+h/K2+h*h/(K1*K2));
co2=dic/(1+K1/h+K1*K2/(h*h));
pH=-log10(h);
pco2=1e6*co2/KH;

%Uncomment to check the calculations
% BOH4=KB*tbor/(h+KB);
% OH=KW/h;
% Ct=(hco3+co3+co2);
% At=(hco3+2*co3+BOH4+OH-h);
end

function [LBC_Data] = LBC_data()
%%
%Val Adrara Ccarb (from van de Schootbrugge et al 2008)
%stratigraphic position, d13C_carb, d18O,
Adrara =  [
           1.8        1.901       -3.056
          2.8         1.87       -2.144
          3.8        1.735       -3.059
          4.8        1.818       -2.774
          5.8        2.138       -2.392
          6.8        2.362       -2.873
          7.8        1.803        -3.07
          8.8        2.059       -1.655
          9.8        2.058       -2.188
           11        1.865       -2.677
         11.6        1.947       -2.812
           12        1.818       -3.548
         12.9        2.205       -2.934
           16        1.943        -3.16
           17        0.672        -3.83
           18        1.643       -3.364
           20        1.358       -2.262
           24        3.276       -4.662
           27         2.54       -2.848
           28        3.195       -2.865
         28.6        3.184       -2.876
         29.3        3.227       -2.874
         30.8         1.41       -3.516
           32        2.236       -2.472
           33        2.068       -2.817
           34        2.016        -2.12
         34.9        1.374       -2.911
         35.9        2.121       -2.908
         36.2        2.388       -2.469
         37.2        2.661       -2.337
         38.3          2.5       -2.799
           40        2.657        -2.88
           44        2.664       -2.644
         44.5        2.635       -2.752
         45.4        2.366       -3.288
         46.2        2.456       -2.909
         46.9        2.544       -2.021
           48        2.418       -2.411
         49.4        2.461       -2.565
         50.3        2.302       -3.176
         52.2        2.193       -3.294
         53.5        2.416       -2.293
           55        2.367       -2.789
         57.2        2.896       -2.289
         57.6        3.145       -1.943
         58.6        2.874       -2.086
         59.2        3.158       -2.268
         59.5         3.23       -3.061
         60.7        4.211       -2.007
         61.2        3.725       -2.967
         61.7        3.817       -2.708
         62.2        3.694       -2.566
         62.7        3.828       -2.435
         63.2        3.746       -2.347
         63.7        3.805       -2.476
         64.2        3.629        -2.53
         64.7        3.786       -2.109
         65.2        3.667       -2.311
         65.7        3.216       -3.041
         66.2         3.27       -3.201
         67.2        3.505       -2.739
         67.7        3.064       -2.943
         68.2        3.303       -2.811
         69.3        3.315       -2.965
         69.8        3.554       -2.444
         70.4         3.11       -2.987
         71.4        3.374       -2.589
         71.8        3.386       -2.985
           79        3.029       -3.398
           80        3.548        -5.52
           81        2.606       -5.664
           82        3.509       -5.243
           83        3.594       -5.449
           84        2.989       -5.592
           85        1.498       -5.839
           86         3.15       -5.835
           87         4.04        -5.38
           88        3.776       -5.353
           89        3.252       -5.862
           90        2.023       -5.893
           91        3.178        -6.05
           92        -0.87        -5.86
           93        3.325       -5.298
           94        3.779       -5.759
           95        4.396       -5.323
           96        4.376       -5.641
           97          4.2       -5.857
           98        4.416        -5.57
        100.6        4.364       -5.681
        101.4        4.781       -5.996
        102.5        4.527       -6.005
        103.5        4.703       -5.577
        105.5        5.043       -4.108
        106.5        5.048       -5.182
          108        5.023       -5.397
          109        5.235       -4.208
          110        5.125       -5.846
          111        5.428       -4.203
          112         5.45       -4.684
        113.4        5.586       -4.477
        114.4        3.005       -5.579
          115        5.722       -4.069
          116        5.464       -4.684
          117        5.544        -4.63
          118        5.558       -4.191
          119        5.492       -4.275
          120        5.656       -3.958
          121        5.231       -4.659
          122        5.456       -4.379
          123        5.559       -3.588
          124        4.536       -3.518
          126        5.759       -3.858
          127          5.5       -4.763
          128        5.213       -4.718
          129        5.435       -4.324
          130        5.503       -4.342
          131        4.249        -4.28
          132        4.378       -5.046
          133        4.784       -4.601
          134         4.98       -4.355
          135        5.212       -3.825
          136        5.412       -3.891
          137        5.398       -4.032
          138        5.392       -4.106
          139        5.209       -4.041
          140        5.374       -4.523
          141        5.174       -4.743
          142         5.46       -3.595
          143        4.803       -4.181
          144        5.206       -4.085
          145        5.067       -4.432
          146        5.118       -4.365
          147        5.189       -4.084
          148        5.039       -3.982
          149        4.223       -5.037
          150        5.066       -4.543
          151        4.814       -4.793
          152        4.599       -5.271
          153        4.832       -3.393
          154        4.838        -3.83
          155        3.573        -4.77
          156        4.578       -3.995
          157        4.732       -3.921
          158        4.946       -3.545
          159        4.704       -4.358
          160        4.902        -3.78
          161        4.424       -4.623
          162        4.271       -5.336
          163        4.771       -4.643
          164        4.682       -5.116
          165        4.577       -4.765
          166        4.878       -3.621
          167        4.444       -4.397
          168        4.112       -4.352
          169        4.438       -4.281
          170        4.649       -3.713
          171        4.598       -4.501
          172        4.654       -4.181
          173        4.695       -3.593
          174        4.658       -3.496
          175        4.756       -3.333
          176        4.654       -3.889
          177        4.603       -4.135
          178        4.662       -3.811
          179        4.606       -3.674
          180        4.058       -3.183
          181        4.146       -4.154
          182        4.351       -4.559
          183        4.387       -3.781
          185         4.44       -3.522
          186        4.327       -4.463
        186.8        4.461       -4.043
          189        4.473       -3.789
          190        4.414       -3.989
          191        4.341       -3.685
          192        4.265       -3.617
          193        4.378       -3.684
          194        4.422       -3.501
          195        4.271       -3.629
          196        4.286       -3.853
          197        4.365       -3.677
          198        4.004       -3.906
          200        4.173       -3.923
          201        4.189       -3.413
          202        4.097         -3.5
          203        4.059       -3.379
          204        4.075       -3.337
          205         4.06        -3.42
          206         4.08       -3.355
          207        4.206       -3.351
          208          4.3        -3.45
          209        4.131       -3.708
          210        4.044       -3.824
          211        4.146       -3.809
          212        3.992       -4.322
          213        3.313       -4.621
          214        4.137       -3.287
          215        4.105       -3.172
          216         4.06       -3.578
          217        3.802       -3.929
          218        3.863       -4.269
          219        4.026        -3.28
          221        3.838       -3.295
          222        3.908       -3.939
          226        3.564       -4.243
          227        3.454       -3.973
          228        3.436        -3.59
          229        3.396       -3.416
          230        3.543       -3.693
          231        3.351       -3.796
          232        2.945       -3.915
          233        3.405       -3.667
          234        3.255       -3.574
          235        3.187       -3.741
          236        3.087       -3.574
          237        3.181       -3.657
          238        3.164       -4.151
          239        2.871       -4.571
          240        3.072       -4.688
        242.1        3.052       -4.257
       242.22        2.481       -4.765
       242.75        3.204         -4.8
       243.03        2.661       -4.115
       243.33        2.378        -4.54
        243.8        2.179       -6.361
        244.4        2.859       -4.447
        244.9        2.519       -4.189
        245.3        2.605       -3.994
          246        2.618       -3.662
        246.5         2.62       -4.366
          247        2.587       -3.741
        247.5        2.522       -3.669
        247.6         2.44       -4.262
          248        2.336       -3.404
        248.5        2.394       -3.697
        249.2        2.409       -3.856
        249.5        2.059       -4.132
        250.2        2.162       -5.036
        250.7        1.779       -4.544
        251.1        2.058       -4.509
        251.4        1.536       -5.192
          252         1.63       -5.354
        252.3        1.451       -4.216
          253        1.558       -3.045
        253.5        1.727       -3.673
          254        1.496       -2.386
          255         1.49        -4.15
        255.2        1.727       -3.211
        255.6        1.615       -3.524
          256        2.009       -4.398
        256.6         1.51       -4.866
          257        1.121       -5.347
        257.5        1.538       -3.978
          258        1.416       -4.049
        258.5        1.792       -4.307
          259        2.107       -3.129
        259.4        1.901       -3.525
          260         1.53       -3.218
        260.5        1.813       -3.326
          261        1.512       -3.209
        261.5        1.539       -5.549
          262         1.87       -3.035
        262.4         1.95       -2.859
          263        1.348       -3.062
        263.5        0.202       -5.715
          264        1.701       -3.489
        264.5        1.588       -3.471
        265.2        1.508       -2.719
        265.6        1.229       -3.913
          266        1.187       -3.951
        266.5        1.095       -3.803
          267        1.281       -4.087
        267.5         1.34       -4.168
        268.1         1.49       -3.611
        268.5        1.433       -4.461
        269.1         1.22       -3.141
        269.5        1.608       -4.798
          270        1.332       -4.314
        270.5        1.313       -3.719
          271        1.161       -3.335
        271.5        1.441       -3.243
          272        0.938       -3.599
        272.8       -0.921       -4.067
          273        1.248       -2.717
        273.5        0.966       -4.084
          274        0.836       -4.848
        274.5        0.786       -3.613
          275        0.727       -4.668
          276        1.031       -2.957
        276.5        1.517       -3.131
          277        0.755       -2.701
        277.5        1.048       -3.638
          278        1.013       -4.167
        278.1        0.821       -4.651
        278.5        0.988       -3.459
          279        1.067       -3.618
        279.5        1.124       -3.642
          280        1.177       -3.793
        280.5        1.108       -3.493
        281.3        2.591       -2.268
        281.5        1.215       -3.353
        282.2        1.239       -3.209
          283        1.162       -2.917
        283.5         1.18       -4.159
          284        0.802       -3.594
        284.6        1.013       -4.115
          285        -0.01       -5.562
        285.8       -0.504       -6.417
        286.5        0.795       -6.027
        287.3       -0.595       -6.876
        287.5        1.115       -2.791
          288       -0.976       -5.357
        288.5        0.063       -4.506
          289        0.376       -4.314
        289.5        0.775       -3.797
        290.5        0.319       -3.028
          291        0.876       -4.175
        291.5        0.856       -4.635
          293        0.989       -4.136
          294        0.969       -4.996
          295       -3.146       -5.411
          296        0.485       -3.628
          297        1.017        -4.62
          299        0.602       -6.434
          302        1.566       -4.535
          305        1.212       -6.217
          308        1.515        -4.34
          311        0.898       -3.974
          314        1.109       -6.684
          317        1.341       -4.335
          320        1.421       -4.041
          323        1.472       -5.782
          324        0.226       -7.752
          326        1.367       -10.76
          327         0.94         -5.5
          329        0.641       -2.957
          332         0.17       -7.549
          335        0.867       -4.169
          338        0.882       -3.804
          344        0.858       -4.905
          347        0.816       -2.928
          349        2.008       -4.049
          350        0.577       -2.968
          353         0.44       -2.966
          356        0.512       -3.752
          359        0.384       -3.163
          362        0.626       -2.701
          365        0.591       -2.964
          368        0.533       -2.634
          371        0.407       -4.145
          374       -0.462       -5.127
          377        0.452       -2.909
          380        0.414       -2.964
          383        0.393       -3.227
          386        0.417       -2.987
          389        0.694        -2.95
          392        0.672       -2.564
          395        0.527       -2.924
          398         0.67       -2.668
          401        0.357       -5.185
          404        0.483       -5.738
          407        0.593       -3.013
          410        0.479       -3.249
          413        0.733       -2.065
          416        0.628       -4.084
          419        0.884       -2.529
          422        0.796       -1.957
        ];
    
   
    
        
% Val Adrara Corg (from Bachan et al 2012):
%stratigraphic elevation, d13C_org, d13C_carb

Adrara_Corg = [
          1.8      -27.078        1.901
          7.8      -28.442        1.803
         11.6      -27.615        1.947
         14.6      -28.191        2.205
           17      -27.802        0.672
           18      -27.445        1.643
           20       -27.81        1.358
           27      -26.255         2.54
         30.8      -26.475         1.41
         46.9      -27.082        2.544
         51.2      -27.673        2.302
           55      -29.602        2.367
         57.2      -28.594        2.896
         62.2      -29.015        3.694
         64.2      -27.294        3.629
         69.3      -29.613        3.315
         71.8      -27.566        3.386
           77      -28.597        3.386
           86      -26.481         3.15
           86        -26.6         3.15
           94      -27.261        3.779
           94      -27.026        3.779
           95      -26.604        4.396
           95      -26.913        4.396
        100.6      -27.024        4.364
          108      -26.857        5.023
          112      -26.955         5.45
          122       -26.52        5.456
          147       -27.21        5.189
          150       -26.86        5.066
          162      -28.066        4.271
          163      -27.202        4.771
          172      -27.524        4.654
          182        -26.8        4.351
          196      -28.168        4.286
          200      -27.921        4.173
          205      -27.453         4.06
          210      -27.641        4.044
          221      -28.125        3.838
          234      -27.606        3.255
          241      -27.964        3.072
        246.6       -29.28         2.62
          251      -29.306        1.779
        255.5      -29.892        1.727
        265.5      -29.301        1.508
        266.5      -29.849        1.095
        271.6      -29.182        1.441
        276.5      -30.297        1.517
        281.2      -30.164        1.108
        284.8      -28.729        1.013
          291      -29.914        0.876
          307      -29.286        1.212
          310       -28.39        1.515
        319.7      -29.204        1.341
        319.7      -29.644        1.341
          325      -28.992        0.226
          331      -30.561        0.641
          343       -29.29        0.882
          356      -27.508        0.512
          364      -29.142        0.626
          370      -29.334        0.533
          376      -26.254       -0.462
          388      -28.239        0.417
          394      -28.401        0.672
          394      -29.142        0.672
          403      -28.608        0.357
          421      -28.713        0.884
          445      -29.671        0.796    
];

%Val Adrara detailed (ADR) Ccarb: 
%stratigraphic position, d13C_carb, d18O,

ADR = [
            0         2.48        -2.74
          0.5         2.59        -2.45
            1          2.6        -2.47
          1.5         2.64         -1.9
            2          2.5        -2.31
          2.5         2.61        -2.61
            3         2.48         -2.7
          3.5          2.4        -2.69
            4         1.87        -2.89
          4.5         2.02        -3.45
            5          2.7        -2.52
          5.5         2.51        -3.02
            6          2.5        -2.64
          6.5         2.57        -2.66
            7         2.37        -2.59
          7.5         2.13         -3.5
            8         2.03        -3.14
            9         1.97           -4
          9.5         2.03        -3.72
           10         1.67        -3.21
         10.5         2.37        -2.91
           11         2.67        -2.47
         11.5         2.24        -3.68
           12         2.22        -2.31
         12.5         1.87        -3.32
           13         2.14        -3.44
           14        -3.23        -5.39
         14.5        -4.55        -5.59
        14.75        -0.73         -5.2
        14.85        -1.48        -5.31
           15         1.82        -3.45
         15.2         2.76        -2.56
         15.5         2.83        -2.28
        15.75         2.99        -2.07
           16         3.37        -1.68
        16.25         3.46        -2.38
        16.75         3.84        -2.53
           17         3.76        -2.13
        17.25         4.02        -2.39
         17.5         3.91        -2.42
        17.75         3.78        -2.58
           18         3.82        -2.33
        18.25         3.82        -2.24
         18.5         3.68        -2.43
        18.75         4.16        -1.72
           19         3.67        -2.62
        19.25         3.73        -2.22
         19.5         3.57        -2.62
        19.75         3.78         -2.2
           20         3.83        -2.11
        20.25          3.6        -2.66
         20.5         3.77        -2.19
        20.75         3.18        -3.35
           21          3.7        -2.09
        21.25         3.57        -2.44
         21.5         3.53        -2.59
        21.75         3.35        -3.31
           22          2.7        -2.47
        22.25         3.55        -2.42
         22.5          3.6        -2.67
        22.75         3.44        -2.51
           23         3.54        -2.69
        23.25         3.49        -2.64
         23.5         3.05        -2.58
        23.75         3.26        -2.66   
];


%Pozzo Glaciale section Ccarb (PZG) 
%         el,       d13C,       d18O

PZG = [   
          0.4         5.55        -6.42
          0.8         5.67        -4.99
          1.9         5.48        -6.48
          2.3         5.46        -5.81
          2.9         5.65        -5.09
          3.5         5.81        -5.74
          4.9         4.97        -6.89
          5.3         5.55        -6.08
            7         5.51        -7.22
          8.9         5.36        -5.13
         11.9         5.33        -5.79
         14.5         6.05        -7.55
         15.1         5.29        -5.15
           16         5.32        -4.55
         16.7         5.27        -6.06
         17.4         5.26        -5.77
         18.1         5.28        -4.74
         18.8         6.04        -8.74
         20.1         4.93         -6.1
           21         5.31        -6.13
         21.7          5.2        -5.65
         22.5         5.12        -6.63
         23.5         5.45        -9.43
         24.2         5.18        -6.04
         25.4         5.15        -5.96
         26.3         5.84        -8.62
           27         4.15        -8.77
         27.5         5.19        -9.25
         28.2            0            0
           29          4.9        -5.04
         29.6         5.03        -4.98
         32.1          5.1        -5.58
           33         4.99        -4.22
           34          5.1        -4.94
         35.5         4.47        -5.71
         37.2         4.97         -5.5
           39         4.99        -5.75
         41.5         4.75        -5.06
         43.7         4.55        -5.12
         46.7         4.62        -4.56
         49.9         4.16        -5.43
         51.5         4.19        -4.93
         54.4            4         -5.2
         57.2         4.33         -5.4
         65.8         3.42        -5.43
           75         2.02        -4.52
         77.9         4.44        -4.74
           79         3.41        -4.57
           81         3.32        -4.93
         83.5         3.28        -4.65
           86         3.16           -7
         91.4         3.07        -5.24
         93.1         2.61        -4.54
         97.3         2.25        -4.86
        100.5         2.16        -4.84
        102.9         1.83        -3.07
        104.6         2.11        -4.18
        106.6         2.03        -3.77
          112         1.89        -3.81
        113.7         0.73        -4.17
        114.6         1.32        -4.37
        115.1         1.34        -3.02
        116.5         1.34        -4.45
        117.1         1.61        -3.96
          118         2.04        -3.84
        120.6         1.74        -3.52
        123.3         1.66        -3.26
        124.8         1.79        -4.46
        126.9         1.78        -4.64
        128.9         2.13        -4.46
        130.8         1.89        -5.29
        131.7         2.08        -3.97
        133.4         2.31        -3.66
        134.2         1.66         -4.8
        135.8         1.64        -4.46
        136.3         1.72        -4.32
        139.3         1.49        -5.03
        140.3         1.36        -4.83
        140.5         0.14        -5.51
        147.9         1.37        -4.83
          150         1.27         -5.5
        155.2        -0.15        -6.11
        153.3         1.45         -5.3
        153.9         0.66        -5.42
        156.9         1.03        -5.86
          159         0.85        -4.11
        160.6         1.37        -4.41
        162.8         1.43        -4.22
        164.2         1.36        -4.18
        165.7         1.38        -5.14
        166.9         1.34        -5.81
          169         1.03        -4.64
        171.3         1.29        -4.39
        175.3        -0.98        -6.45
        177.6         1.72        -4.81
        178.3         1.84        -5.09
        181.5         0.84        -5.89
        186.4          1.8        -5.35
          191         1.55        -6.17
        197.8         1.55        -6.51
          203         0.04        -5.71
        205.4         1.54        -5.55
        208.9         1.56        -4.22
        224.2         1.12        -4.98
        229.5         1.12        -5.06
        234.3         1.56        -4.23
        239.6         1.55        -5.16
        234.6         1.65        -5.79
        245.3         1.02        -5.75
        247.8         0.95        -5.82
          251         1.22        -5.75
        258.8         1.19        -5.05
        264.9          1.2        -4.96
        267.3         1.57        -4.66
        272.5         1.47        -5.09
        276.7         1.75        -3.93
        281.3         0.27        -6.84
        285.4         1.56        -4.45
        292.8         1.33        -5.21
        296.3         1.34        -4.85
        300.5         1.13        -5.53
        309.3          0.7        -5.76
        306.4         0.51        -4.95
        312.6          1.2        -3.91
        315.3         0.74        -5.06
        318.9         0.79        -4.36
        320.4          1.3        -3.19
        323.1         0.88        -3.69
        327.6         0.98        -2.54
        329.6         1.01        -5.21
        336.5         1.02        -3.78
        340.8         0.88        -3.65
          347         0.64        -3.63
        351.5         0.67        -2.96
        354.5         0.82        -2.86
          359         0.73        -2.93
          365         0.94        -2.68
        372.5         0.93        -2.16
        377.5         0.93        -2.58
          382         0.94        -1.72
          388            1        -2.51
          393         0.76        -2.54
          400         1.14        -2.67
        403.4         1.19        -2.44
          406         1.11        -2.26   
];


%Pozzo Glaciale Corg: 
%         el,     d13Corg,     d13Ccarb,  composite el

PZG_Corg = [
          0.4      -26.469       5.5513       126.16
          1.9      -26.586       5.4817       127.89
          2.3      -26.385       5.4624       128.34
          3.4        -25.6       5.6486       129.61
          5.1      -26.735       4.9714       131.56
          5.3      -27.018       5.5463        131.8
         13.4      -26.284       5.3267       141.11
         14.2      -27.475       5.3267       142.03
         15.1      -26.654       5.2866       143.06
         17.3      -26.583       5.2668       145.59
         17.4      -27.057       5.2615       145.71
         18.1      -26.793       5.2814       146.52
         18.8      -25.972        6.041       147.32
         20.1      -25.578       4.9312       148.81
         23.5       -26.78       5.4493       152.72
         25.9       -26.84       5.1507       155.49
         26.3      -26.884       5.8442       155.94
           27      -26.782       4.1523       156.75
           33      -29.056       4.9902       163.65
         36.2      -27.134       4.4665       167.33
           38      -28.553       4.9668        169.4
         40.5      -26.667       4.9872       172.28
         40.5      -26.208       4.9872       172.28
         44.5      -26.704       4.5494       176.88
         49.9      -25.811       4.1589       183.08
         51.1      -27.268       4.1589       184.47
         51.5      -26.883       4.1904       184.93
         57.2        -27.1       4.3347       191.48
           67      -27.122       3.4238       202.75
         75.5      -27.664       2.0238       212.52
         78.1      -28.092       4.4391       215.51
         84.6      -28.621       3.2787       222.99
         87.8       -27.96       3.1628       226.67
         90.3       -27.94       3.1628       229.54
         94.6       -28.62       2.6141       234.49
          112      -30.203       1.8926        254.5
        112.5      -29.369       1.8926       255.07
        117.5      -29.445       1.6058       260.82
        118.7      -28.413       2.0386        262.2
        123.3      -25.863       1.6635        267.5
          130      -29.662       2.1251        275.2
        133.9      -29.846       2.3106       279.69
        135.5      -29.293       1.6633       281.52
        140.3      -29.765       1.3614       287.05
        142.1      -29.135      0.13627       289.12
        153.3      -29.652       1.2734          302
        156.9      -30.011       1.0314       306.13
        159.8      -29.998      0.84759       309.47
        160.6      -30.121       1.3731       310.39
        173.8       -29.61       1.2865       325.57
        177.6      -25.832       1.7242       329.94
        177.6      -28.736       1.7242       329.94
        181.5      -28.598      0.84264       334.43
       183.91      -29.309      0.84264        337.2
        184.9        -25.9      0.84264       338.33
        184.9      -28.317      0.84264       338.33
        190.4      -30.187       1.8013       344.66
        197.8      -29.618       1.5478       353.17
       203.07      -30.019     0.038673       359.23
        205.4      -25.165        1.536       361.91
        212.6      -28.912       1.5596       370.19
        225.7      -29.778       1.1183       385.25
        233.5      -29.725       1.1184       394.22
        249.6      -28.374      0.94674       412.74
          302      -29.695       1.1312          473
        321.8      -29.084       1.3006       495.77
          350      -27.735       0.6355        528.2   
];




%Italcementi Active C_Carb: 
%         el,      d13Ccarb,    d18O
%*updated with "missing" data 48.0- 56.3m Feb 2014
ICA = [
 
          1.5       2.6874      -2.4833
          2.5       2.4007      -2.0539
          3.5       2.4781      -2.2421
         4.05       2.4434      -4.5209
          4.5       2.4693      -3.9683
          4.8       2.8259      -1.9825
            5       2.7966      -1.7982
            6       2.5242      -1.9552
          6.7       2.4685      -1.8304
          7.4       2.4419      -1.3993
          7.5       2.4202      -1.5021
          7.6       2.4201      -1.0665
          7.7       2.3709      -1.5566
          8.3       2.4242      -1.6281
          8.3       2.4293      -1.5643
          8.4       2.3161      -1.5617
          8.5       2.3341      -1.1928
          8.6       2.4285      -1.5474
            9       2.3509      -2.0733
          9.2       2.3849      -1.7267
          9.6       2.8509      -1.6028
          9.7       2.4003      -2.0037
          9.9       2.4752      -1.7679
         10.4       2.4543       -2.039
         10.9       2.4458       -1.955
         11.5       2.4487      -1.6947
         11.6       2.4548      -1.7424
         11.7       2.4007       -1.816
         11.8       2.4542      -1.8542
           12       2.4258      -1.7688
         12.5       3.1837       -2.006
         12.6       2.3692      -2.0441
         12.7       2.4292      -1.7223
         13.2       2.4047      -1.7799
         13.5       2.3835      -1.6151
         13.8       2.3175       -1.508
         14.1       2.2279      -1.5279
         14.3       1.9695      -1.3565
         14.5       2.3262      -1.5149
         14.7       2.2938      -1.1726
           15       2.2592      -1.7122
         15.3       1.9115      -2.2046
         15.6       2.5257      -1.5586
         15.9        2.521      -1.2104
         16.1       2.4783      -1.7926
         16.6       2.4442      -1.7484
         16.8       2.5484      -1.5228
         17.3       2.5667      -1.8639
         17.8       2.6935      -2.0132
         18.7       2.6925      -1.6347
         18.9       2.7482       -1.604
         19.4       2.7584      -1.7941
         20.5       3.0542      -2.5626
         20.6        3.038      -2.3386
         20.7       3.1457      -2.6772
         21.2       3.2588      -1.9363
         21.6        3.307      -1.8539
         21.8       3.3061      -2.0953
           22       3.2084      -2.3142
         22.5       2.9006      -2.2464
         22.9       3.1495      -2.6409
         23.2        3.019      -2.7679
         23.3       2.5603      -1.6885
         23.6        2.911      -2.2995
         23.7        2.563      -1.8222
         23.8       2.5582      -1.7571
         23.9        2.498      -1.5428
         25.2       2.4702      -1.5758
         25.8       2.1522      -1.4256
         26.1       2.3166       -1.688
         26.3       1.3948       -2.116
         26.5       2.2087      -1.4489
         26.9       2.8029      -1.9616
         27.1       1.4149      -1.5076
         27.5       1.4871      -2.8452
         27.7       2.3878       -2.054
         27.9       2.5384      -1.7429
         28.1       2.5596       -2.664
         28.7       2.7292      -1.6759
           29       2.8254      -1.8284
         29.2       2.9446      -1.7898
         29.5        2.782       -1.963
         29.6       2.7333      -1.4628
           30       2.6705      -1.8426
         30.1       2.8567      -1.4672
         30.7       2.8838      -1.3409
         30.9       2.5677      -1.1452
         31.3       2.5334       -1.215
         31.5        2.461      -1.2411
         31.8       2.3941     -0.76453
         32.2       2.3498      -1.2671
         32.6       2.3209      -1.6975
         32.8       2.5389      -1.9612
         33.2       2.5811       -1.986
         33.7        2.677      -2.1886
         33.8       2.7127      -1.6173
         34.1       2.7298      -1.5393
         34.4       2.6511       -1.759
         34.7       2.2831      -1.5046
         35.4       2.5061      -1.5877
         35.6       2.5945      -1.4712
         35.9       2.3223      -1.6888
         36.1       2.5368      -1.9655
         36.3       2.4266      -1.5736
         36.5       2.6846      -1.0861
         36.7       2.7511      -1.8621
         36.9       2.6716      -1.9471
         37.3       2.7524      -1.7456
         37.4       2.7872      -1.2887
         37.8       2.6544      -2.0309
         38.5       2.7033      -2.2892
         38.7       2.7385       -1.865
           39       2.7181      -2.1356
         39.3       3.4216      -1.2366
         39.5       2.8688      -1.9618
         40.1       2.7719      -4.5691
         40.3       2.7179      -2.8594
         40.3       2.6949       -3.863
         40.7       2.6741      -2.6369
         40.7       2.6534      -3.5835
         41.1       2.6848      -4.4815
         41.3       2.6015      -3.8333
         41.5       2.5965      -2.9526
         41.7       2.4476       -3.554
           42       2.3789      -6.3825
         42.3       2.3887      -6.3122
         42.5       2.3386      -6.0528
         42.6       2.3281      -5.4792
         42.8       2.3625       -4.322
         43.1       2.4181      -2.9108
         43.5       2.3569        -2.23
         43.6       2.3834      -2.6714
         43.7       2.3508      -3.4616
         43.8       2.4032      -2.4534
         43.9       2.4714      -2.4459
         44.2       2.2999      -2.5259
         44.5       2.4241      -2.3655
         44.7       2.4819      -2.0149
         44.7       2.5551      -2.0072
         44.7       2.5667      -2.0017
         45.3       2.5769      -3.0771
         45.5       2.5754      -1.9891
           46       2.6699      -2.2155
         46.1       3.1139      -1.8173
         46.5       3.1685      -1.9927
         47.1       3.4888      -1.8823
         47.6       3.9134      -1.7475
         47.8       3.8879      -1.7589
         47.9       3.8965      -1.7728      
             48        3.782       -1.685
         48.2        3.707       -1.815
         48.4        3.971       -1.765
         48.9        4.013       -1.691
         49.2        3.923       -1.818
         49.5        4.096       -1.756
         49.6        3.924       -1.845
         50.2        4.052       -1.799
         51.5        4.343       -1.473
         51.9        3.958       -2.097
         52.2        3.995       -2.149
         52.2        4.033       -2.223
         52.5        4.018        -2.34
         52.7        3.993       -2.175
           53        4.133       -1.843
         53.1        4.004       -2.378
         53.6        3.973       -2.274
         53.9        3.995       -2.477
         54.2        3.989       -2.783
         54.5        3.954       -1.954
         54.8        3.888       -2.192
         55.4        3.945       -2.267
         55.7        3.939       -2.109
           56        3.889        -2.12
            56.3        3.932       -2.066
         56.3       3.7558      -2.0962
         56.6       3.6798      -1.9199
         56.8       3.6675      -1.9943
         57.2        3.624      -1.9856
         57.8       3.3837      -2.1051
         58.3       3.3797       -1.924
         58.5       3.1562      -2.6942
         58.6       3.4853       -2.016
         58.9       3.5426      -1.9967
         59.2       3.4088      -2.0861
         60.7       3.7557      -2.0195
           61       3.7332       -2.188
         61.5       3.9878      -1.9649
           62       3.7745      -2.0711
         62.5       3.7963      -1.9939
         63.5       3.8684      -2.3494
           64       3.7866      -2.0368
         64.5       3.4868      -1.9031
         64.6       3.6512      -2.0043
         65.1       3.6631      -2.0664
         65.3       3.7929      -2.5618
           66       3.4357      -2.5901
         67.4       3.4516      -2.1248
         68.6       3.3826      -2.2063
         69.2       3.2436      -2.4109
         69.7       2.9935      -2.4266
         69.8       3.3436      -1.9581
           70       3.1037      -2.1694
         70.4       2.8021      -2.0679
           71       2.5703      -2.7785
         71.2       2.7694      -2.7655
         71.7       2.6692      -2.2661
         72.1       2.4584      -2.6346
         72.2        2.291      -2.4329
         72.6       3.0471      -2.5383
         72.9       3.3491      -2.1398
           73       3.0995      -2.6691
         73.5       3.0485      -2.6632
         73.8       3.1797      -3.2056
           74       3.4376      -2.7034
         74.8        3.288      -2.4718
         75.2       3.1555      -2.8104
         75.9       3.3426       -2.165
         76.5       3.5445      -2.1497
           77       3.6648      -2.0547
         77.3       3.3259      -2.8653
         77.5       3.5506      -2.1373
           78        3.195      -3.0528
         78.3         3.39      -4.9407
         78.5       3.2436       -2.991
         78.8       3.2357      -3.1736
           79       3.2363      -2.9057
         79.2        3.759      -1.5224
           80       3.3595      -2.5258
         80.5       3.3975      -3.1554
           81        3.467      -3.1833
         81.5        3.287      -3.4505
         82.4       3.5241      -3.6037
         82.5       3.7214       -3.159
         83.2       3.4599      -3.6541
         83.6     -0.12925      -5.4957
         84.2       3.8154      -2.9424
         84.4       1.7686      -4.3907
         84.8       3.8015      -3.0791
         85.5       3.8239      -3.0533
           86       3.8811      -2.7958
         86.5       3.8296       -2.584
         86.9       3.8081      -4.2245
         87.4       3.9424      -2.8481
         87.9       4.0081      -2.6712
         88.5       3.5588      -3.7159
         89.1       3.8018      -5.0786
         90.2       3.6849       -4.459
         90.7       3.8858      -3.5891
         91.6       3.8883      -4.2048
           92       3.8118      -3.8406
         92.4       3.8726      -3.7788
           93       3.9298      -3.9733
         93.5       3.8583      -6.5611
         93.9       3.8955      -4.3658
         95.1       3.8641      -6.2005
         95.5       3.9689      -5.6336
         95.9        3.151      -7.0426
         96.4       3.6377      -5.6487
         97.1       4.1493      -6.1347
         97.6       3.8799      -7.1805
         97.9       1.4499      -6.6833
         98.3       3.8821      -6.6984
         99.2       4.2499      -7.4055
         99.4       3.5424      -7.8108
         99.9        3.779       -8.254
        102.1        3.773      -8.5092
        103.4       4.1539      -8.2081
        104.1       4.1861      -7.9721
        105.1       3.8806      -7.4005   
];

%Italcementi Inactive quarry C_Carb: 
%         el,      d13Ccarb,      d18O

ITX = [
            0         2.57         -6.1
          0.5         2.47        -6.11
          1.1         2.61        -3.18
          1.2         2.62        -2.12
          1.6         2.96        -2.04
          1.7         2.91        -2.16
          1.7         2.82        -2.17
         1.72         3.12        -1.87
         1.73         3.08        -1.92
         1.74         3.28        -1.98
         1.75         3.42        -2.18
            2          3.4        -1.84
          2.3         3.41        -1.74
          2.8         3.51        -1.38
          3.4         3.88        -1.74
          5.3         2.86         -2.8
            6         1.97        -6.82
          6.4         3.09        -5.87
          6.8         2.49        -5.87
          7.4         1.35        -8.72
          7.6         3.76        -2.03
          8.1         3.76         -1.8
          8.7         3.61         -1.8
          9.1         3.58        -2.19
          9.3         3.49        -2.02
          9.4         3.72        -1.92
         10.2         3.58        -2.02
         10.7         3.49        -1.97
         10.9         3.59        -2.07
         11.5         3.49        -2.15
         12.7         3.28        -1.99
         13.6         2.56         -2.1
         14.1         2.24        -2.15
         14.8         2.99        -1.96
         15.4         3.03        -2.61
         15.8         3.12        -3.14
         16.2         3.14        -2.83
         16.2         3.13        -2.68
         16.7         3.38        -2.09
         17.6         3.47        -2.51
         17.7         3.62        -3.66
           18         3.69        -2.46
         18.5         3.34        -2.64
         19.4         3.53         -2.4
           20        -0.53        -6.45
         20.6         3.51         -3.5
         21.2         2.37        -4.14
         21.8         3.07         -3.7
         22.2         3.34        -3.66
           23         3.67        -2.98
         23.3         3.67        -2.48
         23.7         1.38         -7.9
         25.3         2.72        -5.05
           26          3.4         -5.2
         27.2         3.51        -4.62
         27.7         2.77         -7.6
         28.3         3.78        -3.97
         28.9         3.65        -5.24
         29.4         3.68         -7.6
         29.7         3.39        -4.52
         29.8         3.91        -6.41
         30.8         2.33        -7.73
         31.8         3.66        -8.51
         32.8         2.73        -8.63
         33.8         2.73        -9.61
         34.8         3.77         -8.1
         35.8         4.23        -8.28
         36.8         4.38        -8.11
         37.8         4.54        -7.78
         38.8         4.01       -10.32
         39.7         3.98        -7.33   
];


%Adrara Boundary Marl (from Bachan et al 2014 submitted): 
%    strat height   d13Ccarb      d18O       d13Corg
ADR_BM = [
            0         1.69        -3.52       -28.23
            0         1.74        -3.48          NaN
          1.50        -3.47         -5.6       -26.65
          1.50        -3.41        -5.58          NaN
          2.00        -3.47        -5.89       -24.51
          2.00        -4.05        -5.78          NaN
          2.25        -0.73        -5.11       -26.11
          2.25         -0.8        -5.06          NaN
          2.50         1.18         -3.4       -28.25
          2.50         1.13        -3.35          NaN
          2.70         2.61        -2.66       -28.88
          2.70         2.82         -2.4          NaN
];

    
%Italcementi FBM section
%    strat height   d13Ccarb      d18O       d13Corg
ICA_FBM = [
        .10         2.65        -1.97       -29.26
        .15         2.37        -2.52       -26.49
        .15          NaN          NaN       -27.98
        .15          NaN          NaN        -27.8
        .22        -0.48        -3.61       -27.15
        .35          NaN          NaN       -27.92
        .35          NaN          NaN       -27.57
        .45         0.14        -3.66       -27.34
        .45         0.47        -4.03       -27.31
        .55          0.5        -4.05       -26.74
        .55          NaN          NaN       -26.63
        .66         0.36        -3.88       -26.81
        .76         0.99        -2.82       -28.77
       1.09         2.46        -2.49       -27.88
       1.16          2.9        -1.97       -28.56
       1.16          NaN          NaN       -28.31   
];

ICA_SBM = [
       0.01         2.32        -2.43       -28.87
      0.015         1.93        -6.15       -28.12
       0.07        -0.87        -8.17       -25.31
       0.07        -4.55        -7.29       -27.29
       0.07        -0.57        -7.57       -26.48
       0.22          NaN          NaN       -26.99
       0.22         -0.8        -5.86       -26.88
       0.45         -0.1        -5.59        -27.2
       0.81          NaN          NaN       -27.97
       0.81         0.26        -5.24       -28.68
       1.01          NaN          NaN       -27.91
       1.01         2.02        -2.93       -28.07
       1.16         2.01        -3.82          NaN
       1.32         2.76        -2.25        -28.1
];


LBC_Data = struct(...
    'Adrara',Adrara,'Adrara_Corg',Adrara_Corg,...
    'ADR',ADR,'PZG',PZG,'PZG_Corg',PZG_Corg,...
    'ICA',ICA,'ITX',ITX,'ADR_BM',ADR_BM,...
    'ICA_FBM',ICA_FBM,'ICA_SBM',ICA_SBM); 


end

function LBC_carb_lines(LBC_Data)

%Adrara - black
Y1 = LBC_Data.Adrara(:,2);
X1 = LBC_Data.Adrara(:,1);
    line(X1,Y1,'Marker','x','LineStyle','none','MarkerFaceColor','k',...
        'MarkerEdgeColor','k');

%ADR  (Adrara detailed) -black
Y1=LBC_Data.ADR(:,2);
X1=LBC_Data.ADR(:,1)*1.1+40.3;
line(X1,Y1,'Marker','x','LineStyle','none','MarkerFaceColor','k',...
       'MarkerEdgeColor','k');

    
%ADR-BM: Adrara boundary marl   
Y1=LBC_Data.ADR_BM(:,2);
X1=LBC_Data.ADR_BM(:,1)*1.1+40.3 + 13;
line(X1,Y1,'Marker','x','LineStyle','none','MarkerFaceColor','k',...
   'MarkerEdgeColor','k');
   
 

%PZG - green
Y1=LBC_Data.PZG(:,2);
X1=LBC_Data.PZG(:,1)*1.15+(243-102*1.15);
line(X1,Y1,'Marker','+','LineStyle','none','MarkerFaceColor','g',...
    'MarkerEdgeColor','g');

%ICA - blue
Y1=LBC_Data.ICA(:,2);
X1=LBC_Data.ICA(:,1).*25/39+26.34;
line(X1,Y1,'Marker','+','LineStyle','none','MarkerFaceColor','b',...
   'MarkerEdgeColor','b');

%ICA-FBM: Italcementi first boundary marl   
Y1=LBC_Data.ICA_FBM(:,2);
X1=LBC_Data.ICA_FBM(:,1) + 55;
line(X1,Y1,'Marker','x','LineStyle','none','MarkerFaceColor','b',...
   'MarkerEdgeColor','b');


%ICA-SBM: Italcementi second boundary marl   
Y1=LBC_Data.ICA_SBM(:,2);
X1=LBC_Data.ICA_SBM(:,1) + 55;
line(X1,Y1,'Marker','x','LineStyle','none','MarkerFaceColor','b',...
   'MarkerEdgeColor','b');


%ITX -red 
Y1=LBC_Data.ITX(:,2);
X1=LBC_Data.ITX(:,1).*28/21.3+53.57;
line(X1,Y1,'Marker','x','LineStyle','none','MarkerFaceColor','r',...
   'MarkerEdgeColor','r');



end
 
function LBC_org_lines(LBC_Data)

%PZG Corg  
X1 = LBC_Data.PZG_Corg(:,4);
Y1 = LBC_Data.PZG_Corg(:,2);
line(X1,Y1,'Marker','o','LineStyle','none','MarkerFaceColor','none',...
'MarkerEdgeColor','r');

   
%ADRARA Corg
X1 = LBC_Data.Adrara_Corg(:,1);
Y1 = LBC_Data.Adrara_Corg(:,2);
line(X1,Y1,'Marker','o','LineStyle','none','MarkerFaceColor','none',...
       'MarkerEdgeColor','b');

   
%ADR-BM Corg
Y1=LBC_Data.ADR_BM(:,4);
X1=LBC_Data.ADR_BM(:,1)*1.1+40.3 + 13;
line(X1,Y1,'Marker','o','LineStyle','none','MarkerFaceColor','none',...
   'MarkerEdgeColor','k');   
   
   


%ICA-FBM Corg
Y1=LBC_Data.ICA_FBM(:,4);
X1=LBC_Data.ICA_FBM(:,1) + 55;
line(X1,Y1,'Marker','o','LineStyle','none','MarkerFaceColor','none',...
   'MarkerEdgeColor','b');


%ICA-SBM Corg
Y1=LBC_Data.ICA_SBM(:,4);
X1=LBC_Data.ICA_SBM(:,1) + 55;
line(X1,Y1,'Marker','o','LineStyle','none','MarkerFaceColor','none',...
   'MarkerEdgeColor','b');


end




