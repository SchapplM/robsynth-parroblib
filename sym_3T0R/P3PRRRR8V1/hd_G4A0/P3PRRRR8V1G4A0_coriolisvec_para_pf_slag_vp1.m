% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V1G4A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:27
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G4A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:24:50
% EndTime: 2020-08-06 17:24:57
% DurationCPUTime: 6.92s
% Computational Cost: add. (37173->296), mult. (97836->631), div. (6885->9), fcn. (124383->40), ass. (0->278)
t913 = sin(qJ(3,2));
t1024 = pkin(2) * t913;
t920 = cos(qJ(2,2));
t914 = sin(qJ(2,2));
t919 = cos(qJ(3,2));
t977 = t914 * t919;
t853 = pkin(2) * t977 - t920 * pkin(5);
t899 = sin(pkin(3));
t901 = cos(pkin(3));
t835 = t901 * t1024 + t853 * t899;
t1040 = 0.1e1 / t835;
t895 = 0.1e1 / t919;
t1010 = t1040 * t895;
t915 = sin(qJ(3,1));
t921 = cos(qJ(3,1));
t950 = rSges(3,1) * t921 - rSges(3,2) * t915;
t1045 = t950 * m(3);
t951 = rSges(3,1) * t919 - rSges(3,2) * t913;
t1044 = t951 * m(3);
t911 = sin(qJ(3,3));
t917 = cos(qJ(3,3));
t952 = rSges(3,1) * t917 - rSges(3,2) * t911;
t1043 = t952 * m(3);
t912 = sin(qJ(2,3));
t978 = t912 * t917;
t964 = t899 * t978;
t987 = t901 * t911;
t918 = cos(qJ(2,3));
t992 = t899 * t918;
t1042 = 0.1e1 / (-pkin(5) * t992 + (t964 + t987) * pkin(2));
t916 = sin(qJ(2,1));
t976 = t916 * t921;
t962 = t899 * t976;
t984 = t901 * t915;
t922 = cos(qJ(2,1));
t988 = t899 * t922;
t1041 = 0.1e1 / (-pkin(5) * t988 + (t962 + t984) * pkin(2));
t872 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t1039 = 0.2e1 * t872;
t925 = m(2) * rSges(2,1);
t1038 = m(3) * rSges(3,3);
t932 = rSges(3,2) ^ 2;
t933 = rSges(3,1) ^ 2;
t865 = (-t932 + t933) * m(3) - Icges(3,1) + Icges(3,2);
t1037 = -t865 / 0.2e1;
t870 = rSges(3,2) * t1038 - Icges(3,6);
t1036 = -t870 / 0.4e1;
t871 = rSges(3,1) * t1038 - Icges(3,5);
t1035 = t871 / 0.4e1;
t858 = t911 * rSges(3,1) + t917 * rSges(3,2);
t1034 = m(3) * (-t899 * t912 * t858 + t952 * t901);
t859 = t913 * rSges(3,1) + t919 * rSges(3,2);
t1033 = m(3) * (-t899 * t914 * t859 + t951 * t901);
t860 = t915 * rSges(3,1) + t921 * rSges(3,2);
t1032 = m(3) * (-t899 * t916 * t860 + t950 * t901);
t1031 = m(3) * t901;
t908 = legFrame(3,2);
t888 = cos(t908);
t928 = xDP(1);
t1001 = t888 * t928;
t893 = 0.1e1 / t917;
t1011 = t1042 * t893;
t1023 = pkin(2) * t917;
t905 = legFrame(3,1);
t882 = cos(t905);
t885 = sin(t908);
t1005 = t882 * t885;
t902 = legFrame(3,3);
t873 = sin(t902);
t879 = cos(t902);
t898 = sin(pkin(6));
t900 = cos(pkin(6));
t840 = t873 * t900 + t879 * t898;
t843 = -t873 * t898 + t879 * t900;
t876 = sin(t905);
t794 = t840 * t1005 + t876 * t843;
t795 = t843 * t1005 - t876 * t840;
t982 = t901 * t918;
t986 = t901 * t912;
t767 = (-t912 * t794 + t795 * t982) * t1023 + pkin(5) * (t918 * t794 + t795 * t986);
t1008 = t876 * t885;
t800 = t843 * t1008 + t840 * t882;
t803 = t840 * t1008 - t843 * t882;
t768 = -(t800 * t982 - t803 * t912) * t1023 - pkin(5) * (t800 * t986 + t918 * t803);
t788 = (-t912 * t840 + t843 * t982) * t1023 + pkin(5) * (t918 * t840 + t843 * t986);
t926 = xDP(3);
t927 = xDP(2);
t936 = 0.1e1 / pkin(2);
t752 = (-t788 * t1001 + t767 * t926 + t768 * t927) * t936 * t1011;
t1030 = pkin(2) * t752;
t1022 = pkin(2) * t919;
t906 = legFrame(2,1);
t883 = cos(t906);
t909 = legFrame(2,2);
t886 = sin(t909);
t1004 = t883 * t886;
t903 = legFrame(2,3);
t874 = sin(t903);
t880 = cos(t903);
t841 = t874 * t900 + t880 * t898;
t844 = -t874 * t898 + t880 * t900;
t877 = sin(t906);
t796 = t841 * t1004 + t877 * t844;
t797 = t844 * t1004 - t877 * t841;
t981 = t901 * t920;
t985 = t901 * t914;
t769 = (-t914 * t796 + t797 * t981) * t1022 + pkin(5) * (t920 * t796 + t797 * t985);
t1007 = t877 * t886;
t801 = t844 * t1007 + t841 * t883;
t804 = t841 * t1007 - t844 * t883;
t770 = -(t801 * t981 - t804 * t914) * t1022 - pkin(5) * (t801 * t985 + t920 * t804);
t789 = (-t914 * t841 + t844 * t981) * t1022 + pkin(5) * (t920 * t841 + t844 * t985);
t889 = cos(t909);
t999 = t889 * t928;
t753 = (t769 * t926 + t770 * t927 - t789 * t999) * t936 * t1010;
t1029 = pkin(2) * t753;
t897 = 0.1e1 / t921;
t1009 = t1041 * t897;
t1021 = pkin(2) * t921;
t907 = legFrame(1,1);
t884 = cos(t907);
t910 = legFrame(1,2);
t887 = sin(t910);
t1003 = t884 * t887;
t904 = legFrame(1,3);
t875 = sin(t904);
t881 = cos(t904);
t842 = t875 * t900 + t881 * t898;
t845 = -t875 * t898 + t881 * t900;
t878 = sin(t907);
t798 = t842 * t1003 + t878 * t845;
t799 = t845 * t1003 - t878 * t842;
t980 = t901 * t922;
t983 = t901 * t916;
t771 = (-t916 * t798 + t799 * t980) * t1021 + pkin(5) * (t922 * t798 + t799 * t983);
t1006 = t878 * t887;
t802 = t845 * t1006 + t842 * t884;
t805 = t842 * t1006 - t845 * t884;
t772 = -(t802 * t980 - t805 * t916) * t1021 - pkin(5) * (t802 * t983 + t922 * t805);
t790 = (-t916 * t842 + t845 * t980) * t1021 + pkin(5) * (t922 * t842 + t845 * t983);
t890 = cos(t910);
t997 = t890 * t928;
t754 = (t771 * t926 + t772 * t927 - t790 * t997) * t936 * t1009;
t1028 = pkin(2) * t754;
t892 = t917 ^ 2;
t1027 = pkin(2) * t892;
t894 = t919 ^ 2;
t1026 = pkin(2) * t894;
t896 = t921 ^ 2;
t1025 = pkin(2) * t896;
t846 = t898 * t918 + t900 * t986;
t849 = -t898 * t986 + t900 * t918;
t812 = t846 * t879 + t873 * t849;
t945 = t873 * t846 - t849 * t879;
t993 = t899 * t917;
t773 = (t812 * t1005 - t945 * t876) * t911 + t795 * t993;
t776 = (-t812 * t1008 - t945 * t882) * t911 - t800 * t993;
t791 = t812 * t911 + t843 * t993;
t852 = pkin(2) * t978 - t918 * pkin(5);
t834 = pkin(2) * t987 + t852 * t899;
t831 = 0.1e1 / t834;
t758 = (-t791 * t1001 + t773 * t926 + t776 * t927) * t893 * t831;
t1020 = pkin(5) * t758;
t847 = t898 * t920 + t900 * t985;
t850 = -t898 * t985 + t900 * t920;
t813 = t847 * t880 + t874 * t850;
t944 = t874 * t847 - t850 * t880;
t991 = t899 * t919;
t774 = (t813 * t1004 - t944 * t877) * t913 + t797 * t991;
t777 = (-t813 * t1007 - t944 * t883) * t913 - t801 * t991;
t792 = t813 * t913 + t844 * t991;
t759 = (t774 * t926 + t777 * t927 - t792 * t999) * t1010;
t1019 = pkin(5) * t759;
t848 = t898 * t922 + t900 * t983;
t851 = -t898 * t983 + t900 * t922;
t814 = t848 * t881 + t875 * t851;
t943 = t875 * t848 - t851 * t881;
t989 = t899 * t921;
t775 = (t814 * t1003 - t943 * t878) * t915 + t799 * t989;
t778 = (-t814 * t1006 - t943 * t884) * t915 - t802 * t989;
t793 = t814 * t915 + t845 * t989;
t854 = pkin(2) * t976 - t922 * pkin(5);
t836 = pkin(2) * t984 + t854 * t899;
t833 = 0.1e1 / t836;
t760 = (t775 * t926 + t778 * t927 - t793 * t997) * t897 * t833;
t1018 = pkin(5) * t760;
t934 = pkin(5) ^ 2;
t935 = pkin(2) ^ 2;
t969 = t753 * t1024;
t1017 = (-pkin(5) * t969 + (t894 * t935 + t934) * t759) * t759;
t1016 = t758 * t1042;
t1015 = t760 * t1041;
t869 = m(2) * rSges(2,2) - t1038;
t1014 = ((t925 + t1043) * t918 - t912 * t869) * t899;
t1013 = ((t925 + t1044) * t920 - t914 * t869) * t899;
t1012 = ((t925 + t1045) * t922 - t916 * t869) * t899;
t1002 = t888 * t893;
t1000 = t889 * t895;
t998 = t890 * t897;
t996 = t899 * t911;
t995 = t899 * t913;
t994 = t899 * t915;
t990 = t899 * t920;
t979 = t901 * t936;
t967 = t911 * t1020;
t743 = t967 - t1030;
t970 = t911 * t1030;
t725 = (((t901 * t752 + t758 * t992) * t1027 - (t970 - t1020) * t964 + t901 * t743) * t1016 + (t752 * t992 + (t892 * t901 - t911 * t964 - t901) * t758) * t1042 * t1030) * t893;
t740 = -pkin(5) * t970 + (t892 * t935 + t934) * t758;
t731 = (t740 * t979 * t1016 + (-t752 * t852 * t996 + t901 * (t752 * t1027 - t967)) * t831 * t752) * t893;
t734 = (-t743 * t1030 + t740 * t758) * t1042;
t749 = t752 ^ 2;
t755 = t758 ^ 2;
t891 = -m(1) - m(2) - m(3);
t971 = 0.2e1 * m(3);
t975 = -t725 * t1014 - t731 * t1034 - t891 * t734 + ((-t755 * t925 - (t755 + t749) * t1043) * t912 - t758 * t918 * (t858 * t752 * t971 + t758 * t869)) * t899 - t749 * t858 * t1031;
t966 = t913 * t1019;
t744 = t966 - t1029;
t963 = t899 * t977;
t726 = (((t901 * t753 + t759 * t990) * t1026 - (t969 - t1019) * t963 + t901 * t744) * t759 + (t753 * t990 + (t894 * t901 - t913 * t963 - t901) * t759) * t1029) * t1010;
t732 = (t979 * t1017 + (-t753 * t853 * t995 + t901 * (t753 * t1026 - t966)) * t753) * t1010;
t735 = (t744 * t1029 - t1017) * t1040;
t750 = t753 ^ 2;
t756 = t759 ^ 2;
t974 = -t726 * t1013 - t732 * t1033 + t891 * t735 + ((-t756 * t925 - (t756 + t750) * t1044) * t914 - t759 * t920 * (t859 * t753 * t971 + t759 * t869)) * t899 - t750 * t859 * t1031;
t965 = t915 * t1018;
t745 = t965 - t1028;
t968 = t915 * t1028;
t727 = (((t901 * t754 + t760 * t988) * t1025 - (t968 - t1018) * t962 + t901 * t745) * t1015 + (t754 * t988 + (t896 * t901 - t915 * t962 - t901) * t760) * t1041 * t1028) * t897;
t742 = -pkin(5) * t968 + (t896 * t935 + t934) * t760;
t733 = (t742 * t979 * t1015 + (-t754 * t854 * t994 + t901 * (t754 * t1025 - t965)) * t833 * t754) * t897;
t736 = (-t745 * t1028 + t742 * t760) * t1041;
t751 = t754 ^ 2;
t757 = t760 ^ 2;
t973 = -t727 * t1012 - t733 * t1032 - t891 * t736 + ((-t757 * t925 - (t757 + t751) * t1045) * t916 - t760 * t922 * (t860 * t754 * t971 + t760 * t869)) * t899 - t751 * t860 * t1031;
t972 = t932 + t933;
t961 = t911 * t865 * t917;
t960 = t913 * t865 * t919;
t959 = t915 * t865 * t921;
t837 = t870 * t917 + t871 * t911;
t929 = 0.2e1 * qJ(3,3);
t949 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - (0.2e1 * rSges(3,3) ^ 2 + t972) * m(3) / 0.2e1 - Icges(3,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,3);
t958 = 0.4e1 * ((t917 * t1035 + t911 * t1036) * t752 + (t961 / 0.2e1 + (t892 - 0.1e1 / 0.2e1) * t872) * t758) * t752 - t734 * t1014 - (cos(t929) * t1037 + t872 * sin(t929) + t949) * t725 - t837 * t731;
t838 = t870 * t919 + t871 * t913;
t930 = 0.2e1 * qJ(3,2);
t957 = 0.4e1 * ((t919 * t1035 + t913 * t1036) * t753 + (t960 / 0.2e1 + (t894 - 0.1e1 / 0.2e1) * t872) * t759) * t753 + t735 * t1013 - (cos(t930) * t1037 + t872 * sin(t930) + t949) * t726 - t838 * t732;
t839 = t870 * t921 + t871 * t915;
t931 = 0.2e1 * qJ(3,1);
t956 = 0.4e1 * ((t921 * t1035 + t915 * t1036) * t754 + (t959 / 0.2e1 + (t896 - 0.1e1 / 0.2e1) * t872) * t760) * t754 - t736 * t1012 - (cos(t931) * t1037 + t872 * sin(t931) + t949) * t727 - t839 * t733;
t867 = -t972 * m(3) - Icges(3,3);
t955 = t734 * t1034 + t837 * t725 + t867 * t731 + t755 * (t892 * t1039 - t872 + t961);
t954 = -t735 * t1033 + t838 * t726 + t867 * t732 + t756 * (t894 * t1039 - t872 + t960);
t953 = t736 * t1032 + t839 * t727 + t867 * t733 + t757 * (t896 * t1039 - t872 + t959);
t948 = t958 * t893;
t947 = t957 * t895;
t946 = t956 * t897;
t942 = pkin(2) * t996 - t852 * t901;
t941 = pkin(2) * t995 - t853 * t901;
t940 = pkin(2) * t994 - t854 * t901;
t939 = t955 * t1011;
t938 = t954 * t1010;
t937 = t953 * t1009;
t857 = pkin(5) * t916 + t922 * t1021;
t856 = pkin(5) * t914 + t920 * t1022;
t855 = pkin(5) * t912 + t918 * t1023;
t811 = -t898 * t857 + t940 * t900;
t810 = -t898 * t856 + t941 * t900;
t809 = -t898 * t855 + t942 * t900;
t808 = t857 * t900 + t940 * t898;
t807 = t856 * t900 + t941 * t898;
t806 = t855 * t900 + t942 * t898;
t787 = -t808 * t875 + t811 * t881;
t786 = -t807 * t874 + t810 * t880;
t785 = -t806 * t873 + t809 * t879;
t781 = -t890 * t836 + (t808 * t881 + t811 * t875) * t887;
t780 = -t889 * t835 + (t807 * t880 + t810 * t874) * t886;
t779 = -t888 * t834 + (t806 * t879 + t809 * t873) * t885;
t1 = [(t956 * t793 * t998 + t973 * ((t940 * t842 + t845 * t857) * t890 + t887 * t836)) * t833 + (t957 * t792 * t1000 + t974 * ((t941 * t841 + t844 * t856) * t889 + t886 * t835)) * t1040 + (t958 * t791 * t1002 + t975 * ((t942 * t840 + t843 * t855) * t888 + t885 * t834)) * t831 + (-t1000 * t1040 * t789 * t954 - t1002 * t1042 * t788 * t955 - t1041 * t790 * t953 * t998) * t936; (-t778 * t946 + t973 * (t781 * t878 - t884 * t787)) * t833 + (-t777 * t947 + t974 * (t780 * t877 - t883 * t786)) * t1040 + (-t776 * t948 + t975 * (t779 * t876 - t882 * t785)) * t831 + (t768 * t939 + t770 * t938 + t772 * t937) * t936; (-t775 * t946 + t973 * (-t781 * t884 - t878 * t787)) * t833 + (-t774 * t947 + t974 * (-t780 * t883 - t877 * t786)) * t1040 + (-t773 * t948 + t975 * (-t779 * t882 - t876 * t785)) * t831 + (t767 * t939 + t769 * t938 + t771 * t937) * t936;];
taucX  = t1;
