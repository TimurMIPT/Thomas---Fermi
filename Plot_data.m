A = readtable('Золото.txt');
x = A(:, 3);
y = A(:, 4);

disp(class(x));
disp(class(y));

% Преобразуем таблицу в обычный массив
x_num = table2array(x);
arrayY = table2array(y);
% Преобразуем5 все элементы массива в числовой формат
y_num = cellfun(@(y) str2double(strrep(y, ',', '.')), arrayY);

% Выводим результат
disp(x_num)
disp(y_num); 

disp(class(x_num));
disp(class(y_num));

plot(x_num, y_num)





