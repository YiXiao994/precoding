clear all;
settings = init_settings();
beam_Centrals = init_beam_central(settings);
users = init_user_positions(settings, beam_Centrals);

figure;
hold on;
Legend_str = cell(1,settings.num_of_Antenna+settings.num_of_Beams);
for n = 1:settings.num_of_Antenna
    theta = 0:0.5:360;
    edge_Coordinates = beam_Centrals(n,:) + [settings.cell_Radius * cosd(theta);settings.cell_Radius * sind(theta)]';
    plot(edge_Coordinates(:,1),edge_Coordinates(:,2),'LineWidth',2);
    Legend_str{n}=['Antenna',num2str(n)];
end
%legend (antenna_Legend_str,'location','best');
for k = 1:settings.num_of_Beams
   x = users.positions(:,2*k-1);
   y = users.positions(:,2*k);
   scatter(x,y,75,'filled');
   Legend_str{k+settings.num_of_Antenna} = ['Beam',num2str(k)];
end
legend (Legend_str,'location','best');

users = cal_distance_to_satellite(users,settings);
users = cal_angle_to_cell_centre(users,settings, beam_Centrals);
channel_Matrix = init_channel_matrix(users , settings);
result = [];
result1 = [];
result2 = [];
result_mmse = [];
result_zf = [];
result2_mmse = [];
result2_zf = [];
result_group_min = [];
result_group_min1 = [];
result_group_min2 = [];
settings.phase_Error_Standard_Deviation = 5;
for num_of_Users = 1:settings.selected_Users_per_Beam
   H = [];
   for k = 1:settings.num_of_Beams
       begin_index = num_of_Users * (num_of_Users - 1) / 2  + 1;
       H = [H,channel_Matrix(1:settings.num_of_Antenna, (settings.users_per_Beam*(k-1)+begin_index):(settings.users_per_Beam*(k-1)+(begin_index + num_of_Users - 1)))];
   end
   result_temp = optimization_changed(H,settings,num_of_Users);
   result_mmse_temp = MMSE(H, settings, num_of_Users);
   result_zf_temp = zero_force(H, settings, num_of_Users);
   result = [result,result_temp];
   result_mmse = [result_mmse, result_mmse_temp];
   result_zf = [result_zf , result_zf_temp];
%    result_group_min_temp = optimization_changed1(H,settings,num_of_Users);
%    result_group_min = [result_group_min, result_group_min_temp];
end

settings.phase_Error_Standard_Deviation = 15;
for num_of_Users = 1:settings.selected_Users_per_Beam
   H = [];
   for k = 1:settings.num_of_Beams
       begin_index = num_of_Users * (num_of_Users - 1) / 2  + 1;
       H = [H,channel_Matrix(1:settings.num_of_Antenna,settings.users_per_Beam*(k-1)+begin_index:settings.users_per_Beam*(k-1)+(begin_index + num_of_Users - 1))];
   end
   result_temp1 = optimization_changed(H,settings,num_of_Users);
   result1 = [result1,result_temp1];
%    result_group_min_temp1 = optimization_changed1(H,settings,num_of_Users);
%    result_group_min1 = [result_group_min1 , result_group_min_temp1];
end

settings.phase_Error_Standard_Deviation = 30;
for num_of_Users = 1:settings.selected_Users_per_Beam
   H = [];
   for k = 1:settings.num_of_Beams
       begin_index = num_of_Users * (num_of_Users - 1) / 2  + 1;
       H = [H,channel_Matrix(1:settings.num_of_Antenna,settings.users_per_Beam*(k-1)+begin_index:settings.users_per_Beam*(k-1)+(begin_index + num_of_Users - 1))];
   end
   result_temp = optimization_changed(H,settings,num_of_Users);
   result2 = [result2,result_temp];
   result_temp_mmse = MMSE(H,settings,num_of_Users);
   result2_mmse = [result2_mmse,result_temp_mmse];
   result_temp_zf = zero_force(H,settings,num_of_Users);
   result2_zf = [result2_zf,result_temp_zf];
%    result_group_min_temp2 = optimization_changed1(H,settings,num_of_Users);
%    result_group_min2 = [result_group_min2 , result_group_min_temp2 ];
end
users_per_Beam = 1:settings.selected_Users_per_Beam;
mean_sum_rate = [];
mean_sum_rate1 = [];
mean_sum_rate2 = [];
mean_sum_rate_group_min = [];
mean_sum_rate_group_min1 = [];
mean_sum_rate_group_min2 = [];
utility = [];
sum_rate = [];
utility1 = [];
sum_rate1 = [];
utility2 = [];
sum_rate2 = [];
sum_rate_mmse = [];
sum_rate_zf = [];
sum_rate_group_min = [];
sum_weighted_SINR = [];
sum_weighted_SINR_mmse = [];
sum_weighted_SINR_zf = [];
sum_weighted_SINR_group_min = [];
jain_index = [];
jain_index1 = [];
jain_index2 = [];
jain_index_mmse = [];
jain_index_zf = [];
jain_index_group_min = [];
time = [];
time1 = [];
time2 = [];
outage_probability = [];
outage_probability_mmse = [];
outage_probability_zf = [];
outage_probability_group_min2 = [];
for i = users_per_Beam
   utility = [utility , result(i).utility];
   sum_rate = [sum_rate , result(i).sum_Rate];
   sum_rate_mmse = [sum_rate_mmse, result_mmse(i).sum_Rate];
   sum_rate_zf = [sum_rate_zf , result_zf(i).sum_Rate];
  % sum_rate_group_min = [sum_rate_group_min, result_group_min(i).sum_Rate];
   sum_weighted_SINR_temp = [result(i).sum_weighted_SINR];
   sum_weighted_SINR_mmse_temp = [result_mmse(i).sum_weighted_SINR];
   sum_weighted_SINR_zf_temp = [result_zf(i).sum_weighted_SINR];
   sum_weighted_SINR = [sum_weighted_SINR, sum_weighted_SINR_temp];
   sum_weighted_SINR_mmse = [sum_weighted_SINR_mmse, sum_weighted_SINR_mmse_temp];
   sum_weighted_SINR_zf = [sum_weighted_SINR_zf, sum_weighted_SINR_zf_temp];
   jain_index = [jain_index , result(i).jain_Index];
   jain_index_mmse = [jain_index_mmse, result_mmse(i).jain_Index];
   jain_index_zf = [jain_index_zf , result_zf(i).jain_Index];
%   jain_index_group_min = [jain_index_group_min, result_group_min(i).jain_Index];
   mean_sum_rate = [mean_sum_rate , mean(result(i).data_Rate)];
%   time = [time , result(i).time];
   utility1 = [utility1 , result1(i).utility];
   sum_rate1 = [sum_rate1 , result1(i).sum_Rate];
   jain_index1 = [jain_index1 , result1(i).jain_Index];
   mean_sum_rate1 = [mean_sum_rate1 , mean(result1(i).data_Rate)];
%   time1 = [time1 , result1(i).time];
   utility2 = [utility2 , result2(i).utility];
   sum_rate2 = [sum_rate2 , result2(i).sum_Rate];
   mean_sum_rate2 = [mean_sum_rate2 , mean(result2(i).data_Rate)];
  % mean_sum_rate_group_min = [mean_sum_rate_group_min , mean(result_group_min(i).data_Rate)];
  % mean_sum_rate_group_min1 = [mean_sum_rate_group_min1 , mean(result_group_min1(i).data_Rate)];
  % mean_sum_rate_group_min2 = [mean_sum_rate_group_min2 , mean(result_group_min2(i).data_Rate)];
   jain_index2 = [jain_index2 , result2(i).jain_Index];
%   time2 = [time2 , result2(i).time];
   if i == 4
       weighted_sum_SINR_iteration = result(i).weighted_sum_SINR_iteration;
   end
   outage_probability_temp = result2(i).outage_Prob;
   outage_probability_temp_mmse = result2_mmse(i).outage_Prob;
   outage_probability_temp_zf = result2_zf(i).outage_Prob;
  % outage_probability_temp_group_min2 = result_group_min2(i).outage_Prob;
   outage_probability = [outage_probability , mean(vec(outage_probability_temp))];
   outage_probability_mmse = [outage_probability_mmse , mean(vec(outage_probability_temp_mmse))];
   outage_probability_zf = [outage_probability_zf , mean(vec(outage_probability_temp_zf))];
   %outage_probability_group_min2 = [outage_probability_group_min2 , mean(vec(outage_probability_temp_group_min2))];
end
figure;
hold on;
plot(users_per_Beam , utility);
plot(users_per_Beam , utility1);
plot(users_per_Beam , utility2);
xlabel('Users per Beam');
ylabel('Sum Utility')
legend(['Phase uncertainty = ',num2str(2),' deg'],['Phase uncertainty = ',num2str(15),' deg'],['Phase uncertainty = ',num2str(30),' deg']);
title('Utility versus users per beam');

figure;
hold on
plot(users_per_Beam , sum_rate * settings.user_Link_Bandwidth / (10^9));
plot(users_per_Beam, sum_rate_mmse * settings.user_Link_Bandwidth / (10^9));
plot(users_per_Beam , sum_rate_zf * settings.user_Link_Bandwidth / (10^9));
%plot(users_per_Beam , sum_rate_group_min * settings.user_Link_Bandwidth / (10^9))
%plot(users_per_Beam , sum_rate1);
%plot(users_per_Beam , sum_rate2);
xlabel('Users per Beam');
ylabel('Sum Rate (Gbits/s)')
%legend(['Phase uncertainty = ',num2str(2),' deg'],['Phase uncertainty = ',num2str(15),' deg'],['Phase uncertainty = ',num2str(30),' deg']);
legend(['Proposed algorithm'], ['MMSE'],['ZF'])%,['Propose algorithm min SINR']);
%result = optimization_outage(channel_Matrix,settings,7);

figure;
hold on
plot(users_per_Beam , sum_weighted_SINR);
plot(users_per_Beam, sum_weighted_SINR_mmse);
plot(users_per_Beam , sum_weighted_SINR_zf);
%plot(users_per_Beam , sum_rate1);
%plot(users_per_Beam , sum_rate2);
xlabel('Users per Beam');
ylabel('Sum weighted SINR')
%legend(['Phase uncertainty = ',num2str(2),' deg'],['Phase uncertainty = ',num2str(15),' deg'],['Phase uncertainty = ',num2str(30),' deg']);
legend(['Proposed algorithm'], ['MMSE'],['ZF']);
%result = optimization_outage(channel_Matrix,settings,7);
figure;
hold on
plot(users_per_Beam , outage_probability);
plot(users_per_Beam, outage_probability_mmse);
plot(users_per_Beam , outage_probability_zf);
%plot(users_per_Beam , outage_probability_group_min2);
%plot(users_per_Beam , sum_rate1);
%plot(users_per_Beam , sum_rate2);
xlabel('Users per Beam');
ylabel('Outage Probability')
%legend(['Phase uncertainty = ',num2str(2),' deg'],['Phase uncertainty = ',num2str(15),' deg'],['Phase uncertainty = ',num2str(30),' deg']);
legend(['Proposed algorithm'], ['MMSE'],['ZF'])%,['Propose algorithm min SINR']);
figure;
hold on;
plot(users_per_Beam , mean_sum_rate * settings.user_Link_Bandwidth / (10^9));
plot(users_per_Beam , mean_sum_rate1 * settings.user_Link_Bandwidth / (10^9));
plot(users_per_Beam , mean_sum_rate2 * settings.user_Link_Bandwidth / (10^9));
% plot(users_per_Beam , mean_sum_rate_group_min * settings.user_Link_Bandwidth / (10^9));
% plot(users_per_Beam , mean_sum_rate_group_min1 * settings.user_Link_Bandwidth / (10^9));
% plot(users_per_Beam , mean_sum_rate_group_min2 * settings.user_Link_Bandwidth / (10^9));
xlabel('Users per Beam');
ylabel('Individual data rate (Gbits/s)')
legend(['Phase uncertainty = ',num2str(5),' deg'],['Phase uncertainty = ',num2str(15),' deg'],['Phase uncertainty = ',num2str(30),' deg'])%,['Phase uncertainty = ',num2str(5),' deg', ' min SINR'],['Phase uncertainty = ',num2str(15),' deg', ' min SINR'],['Phase uncertainty = ',num2str(30),' deg', ' min SINR']);

legend();
figure;
hold on;
plot(users_per_Beam , jain_index);
plot(users_per_Beam, jain_index_mmse);
plot(users_per_Beam, jain_index_zf)
%plot(users_per_Beam, jain_index_group_min)
%plot(users_per_Beam , jain_index1);
%plot(users_per_Beam , jain_index2);

xlabel('Users per Beam');
ylabel('Jain''s Index')
%legend(['Phase uncertainty = ',num2str(5),' deg'],['Phase uncertainty = ',num2str(15),' deg'],['Phase uncertainty = ',num2str(30),' deg']);
legend(['Proposed algorithm'], ['MMSE'],['ZF'])%,['Propose algorithm min SINR'])
count = 0:10;
figure;
hold on;
plot(count, weighted_sum_SINR_iteration)
xlabel('Iteration times')
ylabel('Weighted sum SINR')
% figure;
% hold on;
% plot(users_per_Beam , time);
% plot(users_per_Beam , time1);
% plot(users_per_Beam , time2);
% xlabel('Users per Beam');
% ylabel('Time consumption (s)');
% legend(['Phase uncertainty = ',num2str(2),' deg'],['Phase uncertainty = ',num2str(10),' deg'],['Phase uncertainty = ',num2str(15),' deg']);
% title('Running time versus users per beam');

% num_of_users1 = 4;
% num_of_users2 = 2;
% bar_result_1 = result(num_of_users1);
% bar_result_2 = result(num_of_users2);
% % figure;
% % hold on;
% % b = bar([bar_result_1.data_Rate',bar_result_2.data_Rate']);
% % label = num2cell(1:settings.num_of_Beams);
% % ch = get(b,'children');
% % set(gca,'XTickLabel',label)
% % xlabel('Beam index')
% % ylabel('Data Rate (bit/Hz)')
% % title('Individual data rate in each user group')
% % legend(['Users per beam = ',num2str(num_of_users1)] , ['Users per beam = ',num2str(num_of_users2)]);
% 
% % figure;
% % hold on;
% % b = bar([bar_result_1.antenna_Power',bar_result_2.antenna_Power']);
% % label = num2cell(1:settings.num_of_Antenna);
% % ch = get(b,'children');
% % set(gca,'XTickLabel',label)
% % xlabel('Beam index')
% % ylabel('Power (Watt)')
% % title('Power distribution of each antenna')
% % legend(['Users per beam = ',num2str(num_of_users1)] , ['Users per beam = ',num2str(num_of_users2)]);
% 
% % figure;
% % hold on;
% % b = bar([bar_result_1.antenna_Power',bar_result_2.antenna_Power']);
% % label = num2cell(1:settings.num_of_Antenna);
% % ch = get(b,'children');
% % set(gca,'XTickLabel',label)
% % xlabel('Beam index')
% % ylabel('Power (Watt)')
% % title('Power distribution of each antenna')
% % legend(['Users per beam = ',num2str(num_of_users1)] , ['Users per beam = ',num2str(num_of_users2)]);
% 
% figure;
% hold on;
% b = bar([mean(bar_result_1.outage_Prob')',mean(bar_result_2.outage_Prob')']);
% label = num2cell(1:settings.num_of_Beams);
% ch = get(b,'children');
% set(gca,'XTickLabel',label)
% xlabel('Beam index')
% ylabel('Mean error probability')
% title('Mean error probability of each user group')
% legend(['Users per beam = ',num2str(num_of_users1)] , ['Users per beam = ',num2str(num_of_users2)]);
% 
% % w_Unoptimized = ones(1 ,settings.num_of_Antenna)
% % W_Unoptimized = repmat(w_Unoptimized/norm(w_Unoptimized),settings.num_of_Beams,1) * sqrt(settings.power_per_Antenna)
% % unoptimized_Result.SINR = [];
% % unoptimized_Result.data_Rate = [];
% % for k = 1:settings.num_of_Beams
% %     w_k = W_Unoptimized(:,k);
% %     SINR_Beam = [];
% %     for q = 1:settings.users_per_Beam
% %        h = channel_Matrix(:,(k-1)*settings.users_per_Beam + q);
% %        signal_Power = square(norm(h' * w_k));
% %        interference_Power = 0;
% %        for l = 1:settings.num_of_Beams
% %          if l ~= k
% %            w_l = W_Unoptimized(:,l);
% %            interference_Power = interference_Power + square(norm(h' * w_l));
% %          end
% %        end   
% %     SINR_User = signal_Power / (interference_Power + settings.noise_Power);
% %     SINR_Beam = [SINR_Beam ; SINR_User];   
% %     end
% %     data_Rate = log2(1 + min(SINR_Beam));
% %     unoptimized_Result.SINR = [unoptimized_Result.SINR,SINR_Beam];
% %     unoptimized_Result.data_Rate = [unoptimized_Result.data_Rate , data_Rate];
% % end
% % 
% % unoptimized_Result.sum_Rate = sum(unoptimized_Result.data_Rate * settings.users_per_Beam);
% % 
% % figure
% % hold on;
