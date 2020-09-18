from django.test import TestCase
from django.test import Client
from django.conf import settings
import os
import shutil
import json
from pals.common import *

class PathwayTestCase(TestCase):

    def test_keypath_get_data(self):
        c = Client()
        with open(settings.TEST_ROOT+'/data/pathway/int_df.csv', 'r') as int_csv:
            with open(settings.TEST_ROOT+'/data/pathway/annotation_df.csv', 'r') as annotation_csv:
                response = c.post('/pals_viewer/pathway/get_data', {
                    'int_csv': int_csv,
                    'annotation_csv': annotation_csv,
                    'data_type': 'upload_file'
                })
                self.assertEqual(response.status_code, 200)

                data = response.json()['data']
                self.assertEqual(len(data['int_df']['columns']), 12)
                self.assertEqual(len(data['annotation_df']['columns']), 1)

    def test_analysis(self):

        # copy test files to target cache directory
        if os.path.exists(settings.MEDIA_ROOT+'/keypath_int_df.csv') == True:
            os.remove(settings.MEDIA_ROOT+'/keypath_int_df.csv')
        if os.path.exists(settings.MEDIA_ROOT+'/keypath_annotation_df.csv') == True:
            os.remove(settings.MEDIA_ROOT+'/keypath_annotation_df.csv')
        shutil.copy2(settings.TEST_ROOT+'/data/pathway/int_df.csv', settings.MEDIA_ROOT+'/keypath_int_df.csv')
        shutil.copy2(settings.TEST_ROOT+'/data/pathway/annotation_df.csv', settings.MEDIA_ROOT+'/keypath_annotation_df.csv')

        # send the request and get response
        c = Client()
        response = c.post('/pals_viewer/pathway/analysis', {
            'keypath_int_df_filename': 'keypath_int_df.csv',
            'keypath_annotation_df_filename': 'keypath_annotation_df.csv',
            'experimental_design': json.dumps({
                "comparisons":[
                    {"case":"beer1","control":"beer2","name":"beer1/beer2"},
                    {"case":"beer3","control":"beer4","name":"beer3/beer4"}
                ],
                "groups":{
                    "beer4":["Beer_4_full3.mzXML","Beer_4_full2.mzXML","Beer_4_full1.mzXML"],
                    "beer3":["Beer_3_full3.mzXML","Beer_3_full2.mzXML","Beer_3_full1.mzXML"],
                    "beer2":["Beer_2_full3.mzXML","Beer_2_full1.mzXML","Beer_2_full2.mzXML"],
                    "beer1":["Beer_1_full2.mzXML","Beer_1_full1.mzXML","Beer_1_full3.mzXML"]
                }
            }),
            'pathway_analysis_method':'PLAGE',
            'database':DATABASE_PIMP_KEGG
        })
        self.assertEqual(response.status_code, 200)

        data = response.json()['data']
        self.assertEqual(data['table']['display_columns'], ["pw_name", "unq_pw_F", "tot_ds_F", "F_coverage", "beer1/beer2 p-value", "beer3/beer4 p-value"])
        self.assertEqual(len(data['table']['rows']), 226)

    def test_show_kegg_diagram(self):
        c = Client()
        response = c.post('/pals_viewer/pathway/show_kegg_diagram', {
            'row': json.dumps({
                "id": "ingenza00008",
                "pw_name":"IBA biosynthesis",
                "beer1/beer2 p-value":0.2736992792,
                "beer3/beer4 p-value":0.1698649217,
                "unq_pw_F":3,
                "tot_ds_F":3,
                "F_coverage":100,
                "undefined":0.1698649217,
                "p-value":0.2736992792,
                "Formula Hits":100
            }),
            'pathway_name': 'IBA biosynthesis',
            'id': 'ingenza00008'
        })
        self.assertEqual(response.status_code, 200)

        data = response.json()['data']
        self.assertEqual(len(data['details']), 5)
        self.assertNotEqual(data['details'][4].find("https://www.genome.jp/kegg/pathway/map/ingenza00008.png"), -1)

    def test_show_reactome_diagram(self):
        c = Client()
        response = c.post('/pals_viewer/pathway/show_reactome_diagram', {
            'row': json.dumps({
                "id":"R-HSA-1483115",
                "pw_name":"Hydrolysis of LPC",
                "beer1/beer2 p-value":0.1925641258,
                "beer3/beer4 p-value":1,
                "unq_pw_F":3,
                "tot_ds_F":1,
                "F_coverage":33.33,
                "undefined":1,
                "p-value":0.1925641258,
                "Formula Hits":33.33
            }),
            'pathway_name': 'Hydrolysis of LPC',
            'id': 'R-HSA-1483115'
        })
        self.assertEqual(response.status_code, 200)

        data = response.json()['data']
        self.assertEqual(len(data['details']), 2)
        self.assertNotEqual(data['details'][2].find("https://reactome.org/ContentService/exporter/diagram/R-HSA-1483115.png"), -1)
    
    def test_gnps_get_data(self):
        c = Client()
        with open(settings.TEST_ROOT+'/data/gnps/metadata_df.csv', 'r') as metadata_csv:
            response = c.post('/pals_viewer/gnps/get_data', {
                'gnps_url': 'https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=0a8432b5891a48d7ad8459ba4a89969f',
                'metadata_csv': metadata_csv,
            })
            self.assertEqual(response.status_code, 200)

            data = response.json()['data']
            self.assertEqual(data['metadata_df']['columns'], ["sample", "group"])
            self.assertEqual(data['metadata_df']['groups'], ["Less than 10", "More than 30"])
    
    def test_gnps_analysis(self):
        # copy test files to target cache directory
        if os.path.exists(settings.MEDIA_ROOT+'/gnps_metadata_df.csv') == True:
            os.remove(settings.MEDIA_ROOT+'/gnps_metadata_df.csv')
        shutil.copy2(settings.TEST_ROOT+'/data/gnps/metadata_df.csv', settings.MEDIA_ROOT+'/gnps_metadata_df.csv')

        # send the request and get response
        c = Client()
        response = c.post('/pals_viewer/gnps/analysis', {
            'gnps_metadata_df_filename': 'gnps_metadata_df.csv',
            'gnps_url': 'https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=0a8432b5891a48d7ad8459ba4a89969f',
            "comparisons": json.dumps([{"name":"test","case":"More than 30","control":"Less than 10"}])
        })
        self.assertEqual(response.status_code, 200)

        data = response.json()['data']
        self.assertEqual(data['table']['display_columns'], ["pw_name", "unq_pw_F", "tot_ds_F", "F_coverage", "test p-value"])
        self.assertEqual(len(data['table']['rows']), 983)
    
    def test_gnps_show_details(self):
        # copy test files to target cache directory
        if os.path.exists(settings.MEDIA_ROOT+'/gnps_load_data.pkl') == True:
            os.remove(settings.MEDIA_ROOT+'/gnps_load_data.pkl')
        shutil.copy2(settings.TEST_ROOT+'/data/gnps/gnps_load_data.pkl', settings.MEDIA_ROOT+'/gnps_load_data.pkl')

        c = Client()
        response = c.post('/pals_viewer/gnps/details', {
            'gnps_load_data_filename': 'gnps_load_data.pkl',
            'row': json.dumps({"id":"2004","pw_name":"Molecular Family #2004","test p-value":1,"unq_pw_F":5,"tot_ds_F":5,"F_coverage":100}),
        })
        self.assertEqual(response.status_code, 200)

        data = response.json()['data']
        self.assertEqual(len(data['details']), 7)
    
    def test_ms2lda_get_data(self):

        c = Client()
        with open(settings.TEST_ROOT+'/data/ms2lda/metadata_df.csv', 'r') as metadata_csv:
            response = c.post('/pals_viewer/ms2lda/get_data', {
                'motif_data_url': 'https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=7c34badae00e43bc87b195a706cf1f43',
                'peak_data_url': 'https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=0a8432b5891a48d7ad8459ba4a89969f',
                'metadata_csv': metadata_csv,
                'data_type': 'url'
            })
            self.assertEqual(response.status_code, 200)

            data = response.json()['data']
            self.assertEqual(data['metadata_df']['columns'], ["sample", "group"])
            self.assertEqual(data['metadata_df']['groups'], ["Less than 10", "More than 30"])
    
    def test_ms2lda_analysis(self):
        # copy test files to target cache directory
        if os.path.exists(settings.MEDIA_ROOT+'/ms2lda_metadata_df.csv') == True:
            os.remove(settings.MEDIA_ROOT+'/ms2lda_metadata_df.csv')
        shutil.copy2(settings.TEST_ROOT+'/data/ms2lda/metadata_df.csv', settings.MEDIA_ROOT+'/ms2lda_metadata_df.csv')

        # send the request and get response
        c = Client()
        response = c.post('/pals_viewer/ms2lda/analysis', {
            'ms2lda_metadata_df_filename': 'ms2lda_metadata_df.csv',
            'peak_data_url': 'https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=0a8432b5891a48d7ad8459ba4a89969f',
            'motif_data_url': 'https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=7c34badae00e43bc87b195a706cf1f43',
            "comparisons": json.dumps([{"name":"test","case":"More than 30","control":"Less than 10"}])
        })
        self.assertEqual(response.status_code, 200)

        data = response.json()['data']
        self.assertEqual(data['table']['display_columns'], ["pw_name", "unq_pw_F", "tot_ds_F", "F_coverage", "test p-value"])
        self.assertEqual(len(data['table']['rows']), 102)
    
    def test_ms2lda_show_details(self):
        # copy test files to target cache directory
        if os.path.exists(settings.MEDIA_ROOT+'/ms2lda_load_data.pkl') == True:
            os.remove(settings.MEDIA_ROOT+'/ms2lda_load_data.pkl')
        shutil.copy2(settings.TEST_ROOT+'/data/ms2lda/ms2lda_load_data.pkl', settings.MEDIA_ROOT+'/ms2lda_load_data.pkl')

        c = Client()
        response = c.post('/pals_viewer/ms2lda/details', {
            'ms2lda_load_data_filename': 'ms2lda_load_data.pkl',
            'row': json.dumps({"id":"motif_194","pw_name":"motif_194","test p-value":1,"unq_pw_F":37,"tot_ds_F":19,"F_coverage":51.35}),
        })
        self.assertEqual(response.status_code, 200)

        data = response.json()['data']
        self.assertEqual(len(data['details']), 7)