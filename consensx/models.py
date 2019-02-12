from django.db import models
import datetime

# Create your models here.
class CSX_upload(models.Model):
    id_code = models.CharField(
        max_length=6, default=None, blank=True, null=True
    )
    PDB_file = models.CharField(
        max_length=40, default=None, blank=True, null=True
    )
    NOE_file = models.CharField(
        max_length=40, default=None, blank=True, null=True
    )
    STR_file = models.CharField(
        max_length=40, default=None, blank=True, null=True
    )
    bme_weights_file = models.CharField(
        max_length=40, default=None, blank=True, null=True
    )
    karplus = models.IntegerField(default=1)
    superimpose = models.BooleanField(default=False)
    timestamp = models.DateTimeField(auto_now_add=True)
    fit_range = models.CharField(
        max_length=10, default=None, blank=True, null=True
    )
    r3average = models.BooleanField(default=False)
    svd_enable = models.BooleanField(default=False)
    rdc_lc = models.CharField(
        max_length=6, default=None, blank=True, null=True
    )

    class Meta:
        verbose_name = 'Calculation'
        verbose_name_plural = 'Calculations'

    def __str__(self):
        my_time = self.timestamp.strftime('%Y-%m-%d %H:%M:%S')
        return my_time + " --- " + self.id_code

class CSX_calculation(models.Model):
    id_code = models.CharField(
        max_length=6, default=None, blank=True, null=True
    )
    html_content = models.TextField(default=None, blank=True, null=True)
    timestamp = models.DateTimeField(auto_now_add=True)


    class Meta:
        verbose_name = 'Calculation result'
        verbose_name_plural = 'Calculation results'

    def __str__(self):
        my_time = self.timestamp.strftime('%Y-%m-%d %H:%M:%S')
        return my_time + " --- " + self.id_code

    def returnHTML(self):
        return self.html_content
