function temp_schedule = compute_cooling_schedule(init_temp, final_temp, anneal_rate)

temp_schedule = [];
i = 1;
temp_schedule(i)=init_temp;
while temp_schedule(i) > final_temp
  i = i + 1;
  temp_schedule(i)=temp_schedule(i-1)*anneal_rate;
end

  


