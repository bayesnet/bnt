function schedule = mk_schedule(init_temp, final_temp, anneal_rate)

init_temp = 10; final_temp = 1e-2; anneal_rate = 0.8;
schedule = [];
temp=init_temp;
schedule = [schedule temp];
while temp > final_temp
  temp = temp * anneal_rate;
  schedule = [schedule temp];
end
