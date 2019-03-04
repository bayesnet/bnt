function data = get_data_from_mysql(sql,datasource,username,password,url)

% this function can connect database and list data,its encapsulating
% Matlab function 'database'.
% ps: this function can't work if data size bigger than 1000000.（my mac
%     pro SDRAM memory is 16G)
%     my database project url is https://github.com/WANGXin1205/uci_database
%     you can use idea open the project and save some data,i support excel
%     in this project which downland form uci
%    （http://archive.ics.uci.edu/ml/datasets.html）.
%
% input:
%     sql: is short for structured query language ,like 'select * from
%          emp;'
%     datasource: your datasource name, just support mysql database,
%                 my datasource is 'growlithe'
%     username: your database username,like 'growlithe'
%     password: your database password,like 'growlithe'
%     url: like url = 'jdbc:mysql://localhost:3306/growlithe';
% output:
%     data: cell data
%
% if you find bug, you can send email to me.
%
% Written by WANGXin(growlithe1205@gmail.com)
%
% 这个函数可以连接数据库并获取数据，它是Matlab自带函数'database'的封装。
% 注意: 这个函数在获取数据量超过1000000时不能正常工作。（我的mac pro内存大小为16G）
%      我的数据库项目地址: https://github.com/WANGXin1205/uci_database
%      你可以使用idea打开这个项目并且保存一些数据，我在这个项目里提供了一些从uci
%     （http://archive.ics.uci.edu/ml/datasets.html）下载的excel。
%
% 输入:
%    sql: sql是结构化查询语言的简称，像'select * from emp;'
%    datasource: 数据源，只支持mysql数据库，我的数据源是'growlithe'
%    username: 数据库名称，像'growlithe'
%    password: 数据库密码，像'growlithe'
%    url: 像 url = 'jdbc:mysql://localhost:3306/growlithe';
% 输出:
%    data: cell 类型的 data
%
% 如果你发现了bug，你可以给我发email
%
% Written by WANGXin(growlithe1205@gmail.com)

check_data_in_get_data_from_mysql(sql,datasource,username,password,url)

% default drive
drive = 'com.mysql.cj.jdbc.Driver';

try
    conn = database(datasource,username,password,drive,url);
    curs = exec(conn, sql);
    curs = fetch(curs);
    data = curs.Data;
catch
    disp('check you datasource、username、password、drive、url or sql is correct')
end

end


function check_data_in_get_data_from_mysql(sql,datasource,username,password,url)

% check all varargin are legal

if isempty(sql)
    error('sql is empty')
end

% check sql is select data, not update, delete, insert
sql_words = split(sql);
if ~isequal('SELECT',upper(sql_words{1}))
    error('the sql need select data')
end

if isempty(datasource)
    error('datasource is empty')
end

if isempty(username)
    error('username is empty')
end

if isempty(password)
    error('password is empty')
end

if isempty(url)
    error('url is empty')
end

end