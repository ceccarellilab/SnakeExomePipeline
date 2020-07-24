from requests import get
from requests.exceptions import RequestException
from contextlib import closing
from bs4 import BeautifulSoup
import csv 
import numpy as np

def simple_get(url):
    """
    Attempts to get the content at `url` by making an HTTP GET request.
    If the content-type of response is some kind of HTML/XML, return the
    text content, otherwise return None.
    """
    try:
        with closing(get(url, stream=True)) as resp:
            if is_good_response(resp):
                return resp.content
            else:
                return None

    except RequestException as e:
        log_error('Error during requests to {0} : {1}'.format(url, str(e)))
        return None


def is_good_response(resp):
    """
    Returns True if the response seems to be HTML, False otherwise.
    """
    content_type = resp.headers['Content-Type'].lower()
    return (resp.status_code == 200 
            and content_type is not None 
            and content_type.find('html') > -1)


def log_error(e):
    """
    It is always a good idea to log errors. 
    This function just prints them, but you can
    make it do anything.
    """
    print(e)


list_id = []
with open ('EntrezEnsebml.txt', 'r') as f:
    for row in csv.reader(f,delimiter='\t'):
            list_id.append(row)

allData = []

with open('cbk.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    #for i in range(30):
    for elem in list_id:
        print(elem[1])
        print(elem[2])
        url = 'https://ckb.jax.org/gene/show?geneId=%s&tabType=GENE_LEVEL_EVIDENCE' % (elem[1])
        print(url)
        raw_html = simple_get(url)

        html = BeautifulSoup(raw_html, 'html.parser')

        tables = html.findAll("table")

        if len(tables)>2:
            tables[2]

            allTables = []
                
            table = tables[2]

            #for table in tables:
            output_rows = []
            for table_row in table.findAll('tr'):
                columns = table_row.findAll('td')
                output_row = [elem[2]]
                for column in columns:
                    output_row.append(column.text.strip())
                output_rows.append(output_row)
                tsv_writer.writerow(output_row)   
            print (output_rows)


    

