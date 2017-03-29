## Makefile.am in applications/ash3d/cloud_tester

bin_PROGRAMS += applications/ash3d/cloud_tester/cloud

applications_ash3d_cloud_tester_cloud_SOURCES = \
	applications/ash3d/cloud_tester/cloud_user.cpp \
	applications/ash3d/cloud_tester/cloud_user.h \
	applications/ash3d/cloud_tester/cloud.cpp 

applications_ash3d_cloud_tester_cloud_CPPFLAGS = \
    $(AM_CPPFLAGS) $(FCLAW_ASH3D_CPPFLAGS)

applications_ash3d_cloud_tester_cloud_LDADD = \
    $(FCLAW_ASH3D_LDADD) $(LDADD)

##
