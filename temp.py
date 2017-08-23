import html2text

with open('boring-feeds.html') as f:
    l = f.readlines()

for i in range(len(l)):
    this_line = html2text.html2text(l[i])
    this_name = this_line[this_line.find('[')+1:this_line.find(']')]
    this_url = this_line[this_line.find('(')+1:this_line.find(')')]
    l[i] = this_name + ' | ' + this_url + '\n'


with open('boring-feeds-2.html','w') as f:
    for i in range(len(l)):
        f.writelines(l[i])


