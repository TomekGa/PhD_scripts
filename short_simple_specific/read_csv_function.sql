-- Linux
-- $ sudo apt-get install postgresql-plpython3-<postgres version number>

--CREATE EXTENSION plpython3u WITH CASCADE;

CREATE FUNCTION read_csv_int(header text, column_number integer, url_link text)
    RETURNS SETOF integer
AS $$
    import urllib.request as ur
    data = ur.urlopen(url_link)
    new_data = data.read().decode('utf-8')
    new_data = new_data.strip()
    new_data = new_data.split('\n')
    newer_data = []
    for i in range(len(new_data)):
        if i == 0 and header == 'True':
            continue
        splited = new_data[i].split(',')[column_number-1]
        splited = int(splited)
        newer_data.append(splited)
    return newer_data
$$ LANGUAGE plpython3u;

CREATE FUNCTION read_csv_txt(header text, column_number integer, url_link text)
    RETURNS SETOF text
AS $$
    import urllib.request as ur
    data = ur.urlopen(url_link)
    new_data = data.read().decode('utf-8')
    new_data = new_data.strip()
    new_data = new_data.split('\n')
    newer_data = []
    for i in range(len(new_data)):
        if i == 0 and header == 'True':
            continue
        splited = new_data[i].split(',')[column_number-1]
        splited = str(splited)
        newer_data.append(splited)
    return newer_data
$$ LANGUAGE plpython3u;
