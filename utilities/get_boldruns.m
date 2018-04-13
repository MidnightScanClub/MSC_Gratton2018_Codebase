function [bold_runs,vclist_runs] = get_boldruns(task,sub_id,sess,vcnum)
% A script that allows me to centrally handle how we access bold run names
% and exceptions


% the standard case
switch task
    case 'motor'
        bold_runs = {'1','2'};
        vclist_runs = {vcnum,vcnum};
    case 'mixed'
        bold_runs = {'1','2'};
        vclist_runs = {vcnum,vcnum};
    case 'mem'
        bold_runs = {'_face','_scene','_word'};
        vclist_runs = {vcnum,vcnum,vcnum};
end
    
% the exceptions
if strcmp(sub_id,'MSC02')
    if (sess == 3) && strcmp(task,'mixed')
        bold_runs = {'2','2'};
        vclist_runs = {vcnum,'vc38538'};
    elseif (sess == 4) && strcmp(task,'mixed')
        bold_runs = {'3','4'};
        vclist_runs = {'vc38538','vc38538'};
    elseif (sess == 4) && strcmp(task,'mem')
        bold_runs = {'_face','_scene2','_word'};
        vclist_runs = {'vc38538','vc38538','vc38538'};
    elseif (sess == 8) && strcmp(task,'mem')
        bold_runs = {'_face','_scene','_word'};
        vclist_runs = {vcnum,'vc38538',vcnum};
    end
% elseif strcmp(sub_id,'MSC07')
%     if (sess == 1) && strcmp(task,'mixed')
%         bold_runs = {'2'};
%         vclist_runs = {vcnum};
%     end
elseif strcmp(sub_id,'MSC09')
    if (sess == 6) && strcmp(task,'motor')
        bold_runs = {'2','3'};
        vclist_runs = {vcnum, vcnum};
    end
elseif strcmp(sub_id,'MSC10')
    if (sess == 6) && strcmp(task,'mem')
        bold_runs = {'_face','_word'};
        vclist_runs = {'vc39619','vc39619'};
    end
end

end