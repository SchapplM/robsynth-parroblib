% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRRR12V1G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x14]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR12V1G2A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:50
% EndTime: 2020-08-06 18:24:55
% DurationCPUTime: 4.21s
% Computational Cost: add. (12977->334), mult. (21903->682), div. (2289->18), fcn. (15576->18), ass. (0->288)
t1060 = 2 * qJ(2,3);
t842 = legFrame(3,2);
t812 = sin(t842);
t815 = cos(t842);
t858 = xDP(2);
t859 = xDP(1);
t786 = t812 * t859 + t815 * t858;
t851 = cos(qJ(3,3));
t1025 = t786 * t851;
t846 = sin(qJ(1,3));
t852 = cos(qJ(1,3));
t857 = xDP(3);
t923 = t812 * t858 - t815 * t859;
t777 = -t923 * t846 + t852 * t857;
t845 = sin(qJ(3,3));
t826 = 0.1e1 / t845;
t1035 = t777 * t826;
t1059 = (-pkin(5) - pkin(6));
t825 = pkin(1) - t1059;
t1036 = t777 * t825;
t1049 = t845 * pkin(3);
t805 = qJ(2,3) + t1049;
t796 = 0.1e1 / t805 ^ 2;
t1037 = t777 * t796;
t1007 = t851 * qJ(2,3);
t838 = t851 ^ 2;
t1074 = -t838 + 0.1e1;
t804 = t857 * t825;
t818 = pkin(3) * t857;
t1054 = pkin(3) * t852;
t956 = (t851 + 0.1e1) * (t851 - 0.1e1) * t1054;
t759 = ((t923 * qJ(2,3) + t804) * t852 + (qJ(2,3) * t857 - t923 * t825) * t846 + pkin(3) * t1025) * t845 - t923 * t956 + t786 * t1007 + t1074 * t846 * t818;
t1043 = t759 * t826;
t795 = 0.1e1 / t805;
t797 = t795 * t796;
t744 = (0.2e1 * t796 * t1025 - t759 * t797) * t1035 - (-t1036 + t1043) * t795 * t1037;
t860 = pkin(1) + pkin(5);
t783 = t786 ^ 2;
t868 = 0.1e1 / pkin(3) ^ 2;
t1028 = t783 * t868;
t878 = t845 ^ 2;
t827 = 0.1e1 / t878;
t986 = t827 * t1028;
t935 = t860 * t986;
t968 = 0.2e1 * t759 * t1035;
t1077 = t744 * t1060 + t796 * t968 + t935;
t1062 = 2 * qJ(2,2);
t843 = legFrame(2,2);
t813 = sin(t843);
t816 = cos(t843);
t788 = t813 * t859 + t816 * t858;
t853 = cos(qJ(3,2));
t1024 = t788 * t853;
t848 = sin(qJ(1,2));
t854 = cos(qJ(1,2));
t922 = t813 * t858 - t816 * t859;
t778 = -t922 * t848 + t854 * t857;
t847 = sin(qJ(3,2));
t830 = 0.1e1 / t847;
t1032 = t778 * t830;
t1033 = t778 * t825;
t1048 = t847 * pkin(3);
t806 = qJ(2,2) + t1048;
t799 = 0.1e1 / t806 ^ 2;
t1034 = t778 * t799;
t1006 = t853 * qJ(2,2);
t839 = t853 ^ 2;
t1073 = -t839 + 0.1e1;
t1052 = pkin(3) * t854;
t955 = (t853 + 0.1e1) * (t853 - 0.1e1) * t1052;
t760 = ((t922 * qJ(2,2) + t804) * t854 + (qJ(2,2) * t857 - t922 * t825) * t848 + pkin(3) * t1024) * t847 - t922 * t955 + t788 * t1006 + t1073 * t848 * t818;
t1042 = t760 * t830;
t798 = 0.1e1 / t806;
t800 = t798 * t799;
t745 = (0.2e1 * t799 * t1024 - t760 * t800) * t1032 - (-t1033 + t1042) * t798 * t1034;
t784 = t788 ^ 2;
t1027 = t784 * t868;
t881 = t847 ^ 2;
t831 = 0.1e1 / t881;
t983 = t831 * t1027;
t934 = t860 * t983;
t967 = 0.2e1 * t760 * t1032;
t1076 = t745 * t1062 + t799 * t967 + t934;
t1064 = 2 * qJ(2,1);
t844 = legFrame(1,2);
t814 = sin(t844);
t817 = cos(t844);
t790 = t814 * t859 + t817 * t858;
t855 = cos(qJ(3,1));
t1023 = t790 * t855;
t850 = sin(qJ(1,1));
t856 = cos(qJ(1,1));
t921 = t814 * t858 - t817 * t859;
t779 = -t921 * t850 + t856 * t857;
t849 = sin(qJ(3,1));
t834 = 0.1e1 / t849;
t1029 = t779 * t834;
t1030 = t779 * t825;
t1047 = t849 * pkin(3);
t807 = qJ(2,1) + t1047;
t802 = 0.1e1 / t807 ^ 2;
t1031 = t779 * t802;
t1005 = t855 * qJ(2,1);
t840 = t855 ^ 2;
t1072 = -t840 + 0.1e1;
t1050 = pkin(3) * t856;
t954 = (t855 + 0.1e1) * (t855 - 0.1e1) * t1050;
t761 = ((t921 * qJ(2,1) + t804) * t856 + (qJ(2,1) * t857 - t921 * t825) * t850 + pkin(3) * t1023) * t849 - t921 * t954 + t790 * t1005 + t1072 * t850 * t818;
t1041 = t761 * t834;
t801 = 0.1e1 / t807;
t803 = t801 * t802;
t746 = (0.2e1 * t802 * t1023 - t761 * t803) * t1029 - (-t1030 + t1041) * t801 * t1031;
t785 = t790 ^ 2;
t1026 = t785 * t868;
t884 = t849 ^ 2;
t835 = 0.1e1 / t884;
t980 = t835 * t1026;
t933 = t860 * t980;
t966 = 0.2e1 * t761 * t1029;
t1075 = t746 * t1064 + t802 * t966 + t933;
t828 = t826 * t827;
t1071 = t783 * (t828 * t838 + t826);
t832 = t830 * t831;
t1070 = t784 * (t832 * t839 + t830);
t836 = t834 * t835;
t1069 = t785 * (t836 * t840 + t834);
t1068 = (0.2e1 * t838 - 0.1e1) * t826;
t1067 = (0.2e1 * t839 - 0.1e1) * t830;
t1066 = (0.2e1 * t840 - 0.1e1) * t834;
t1010 = t826 * t851;
t1021 = t795 * t851;
t1061 = -2 * qJ(2,3);
t861 = qJ(2,3) ^ 2;
t866 = pkin(3) ^ 2;
t870 = pkin(1) ^ 2;
t917 = (2 * t1059 * pkin(1)) - (pkin(6) ^ 2) - t866 - t870 + ((-2 * pkin(6) - pkin(5)) * pkin(5));
t989 = t825 * t1043;
t867 = 0.1e1 / pkin(3);
t993 = qJ(2,3) * t786 * t867;
t889 = -t777 * ((-t786 * t825 * t1010 + (t989 + (t1049 * t1061 + t838 * t866 - t861 + t917) * t777) * t795) * t796 + t797 * t989) + ((t1021 * t1036 + t786 * t826) * t845 + t826 * t993) * t827 * t795 * t786;
t1009 = t830 * t853;
t1019 = t798 * t853;
t1063 = -2 * qJ(2,2);
t862 = qJ(2,2) ^ 2;
t988 = t825 * t1042;
t994 = qJ(2,2) * t788 * t867;
t888 = -t778 * ((-t788 * t825 * t1009 + (t988 + (t1048 * t1063 + t839 * t866 - t862 + t917) * t778) * t798) * t799 + t800 * t988) + ((t1019 * t1033 + t788 * t830) * t847 + t830 * t994) * t831 * t798 * t788;
t1008 = t834 * t855;
t1017 = t801 * t855;
t1065 = -2 * qJ(2,1);
t863 = qJ(2,1) ^ 2;
t987 = t825 * t1041;
t995 = qJ(2,1) * t790 * t867;
t887 = -t779 * ((-t790 * t825 * t1008 + (t987 + (t1047 * t1065 + t840 * t866 - t863 + t917) * t779) * t801) * t802 + t803 * t987) + ((t1017 * t1030 + t790 * t834) * t849 + t834 * t995) * t835 * t801 * t790;
t1058 = pkin(1) * t744;
t1057 = pkin(1) * t745;
t1056 = pkin(1) * t746;
t1055 = pkin(3) * t851;
t1053 = pkin(3) * t853;
t1051 = pkin(3) * t855;
t1046 = t744 * t826;
t1045 = t745 * t830;
t1044 = t746 * t834;
t774 = t777 ^ 2;
t1040 = t774 * t796;
t775 = t778 ^ 2;
t1039 = t775 * t799;
t776 = t779 ^ 2;
t1038 = t776 * t802;
t732 = t889 - t1058;
t1022 = t795 * t732;
t734 = t888 - t1057;
t1020 = t798 * t734;
t736 = t887 - t1056;
t1018 = t801 * t736;
t1016 = t812 * t846;
t1015 = t813 * t848;
t1014 = t814 * t850;
t1013 = t815 * t846;
t1012 = t816 * t848;
t1011 = t817 * t850;
t926 = t777 * t795 * t993;
t984 = t838 * t1028;
t944 = t828 * t984;
t1004 = -0.2e1 * t926 * t1010 + t1077 * t845 + t860 * t944;
t1003 = 0.2e1 * t926 + (-t935 + t1077) * t851;
t925 = t778 * t798 * t994;
t981 = t839 * t1027;
t943 = t832 * t981;
t1002 = -0.2e1 * t925 * t1009 + t1076 * t847 + t860 * t943;
t1001 = 0.2e1 * t925 + (-t934 + t1076) * t853;
t924 = t779 * t801 * t995;
t978 = t840 * t1026;
t942 = t836 * t978;
t1000 = -0.2e1 * t924 * t1008 + t1075 * t849 + t860 * t942;
t999 = 0.2e1 * t924 + (-t933 + t1075) * t855;
t992 = t744 * t1010;
t991 = t745 * t1009;
t990 = t746 * t1008;
t829 = 0.1e1 / t878 ^ 2;
t985 = t783 * t829 * t851;
t833 = 0.1e1 / t881 ^ 2;
t982 = t784 * t833 * t853;
t837 = 0.1e1 / t884 ^ 2;
t979 = t785 * t837 * t855;
t977 = t795 * t1016;
t976 = t795 * t1013;
t975 = t798 * t1015;
t974 = t798 * t1012;
t973 = t801 * t1014;
t972 = t801 * t1011;
t971 = t852 * t795 * t744;
t970 = t854 * t798 * t745;
t969 = t856 * t801 * t746;
t965 = -0.2e1 * t777 * t1016;
t964 = 0.2e1 * t777 * t1013;
t963 = -0.2e1 * t778 * t1015;
t962 = 0.2e1 * t778 * t1012;
t961 = -0.2e1 * t779 * t1014;
t960 = 0.2e1 * t779 * t1011;
t953 = t744 * t845 * t1021;
t952 = t745 * t847 * t1019;
t951 = t746 * t849 * t1017;
t947 = t786 * t852 * t1037;
t946 = t788 * t854 * t1034;
t945 = t790 * t856 * t1031;
t941 = t744 * t977;
t940 = t744 * t976;
t939 = t745 * t975;
t938 = t745 * t974;
t937 = t746 * t973;
t936 = t746 * t972;
t932 = qJ(2,1) * t1038 + t860 * t746 - t887;
t931 = qJ(2,2) * t1039 + t860 * t745 - t888;
t930 = qJ(2,3) * t1040 + t860 * t744 - t889;
t929 = t846 * t953;
t928 = t848 * t952;
t927 = t850 * t951;
t768 = -t986 - t1040;
t920 = -t829 * t984 + t768;
t769 = -t983 - t1039;
t919 = -t833 * t981 + t769;
t770 = -t980 - t1038;
t918 = -t837 * t978 + t770;
t916 = t795 * t1071;
t915 = t798 * t1070;
t914 = t801 * t1069;
t913 = t930 * t1010;
t912 = t931 * t1009;
t911 = t932 * t1008;
t907 = (t828 * t1028 + t768 * t826) * t851;
t906 = (t832 * t1027 + t769 * t830) * t853;
t905 = (t836 * t1026 + t770 * t834) * t855;
t904 = t796 * (-t774 * t812 + t786 * t964);
t903 = t796 * (-t774 * t815 + t786 * t965);
t902 = t799 * (-t775 * t813 + t788 * t962);
t901 = t799 * (-t775 * t816 + t788 * t963);
t900 = t802 * (-t776 * t814 + t790 * t960);
t899 = t802 * (-t776 * t817 + t790 * t961);
t793 = qJ(2,3) * t852 - t825 * t846;
t766 = (t812 * t1055 - t793 * t815) * t845 + t815 * t956 + t812 * t1007;
t898 = (t759 * t964 - t766 * t774) * t797;
t772 = (t815 * t1055 + t793 * t812) * t845 + t1074 * t812 * t1054 + t815 * t1007;
t897 = (t759 * t965 - t772 * t774) * t797;
t780 = t846 * t805 + t825 * t852;
t896 = (-t774 * t780 + t852 * t968) * t797;
t794 = qJ(2,2) * t854 - t825 * t848;
t767 = (t813 * t1053 - t794 * t816) * t847 + t816 * t955 + t813 * t1006;
t895 = (t760 * t962 - t767 * t775) * t800;
t773 = (t816 * t1053 + t794 * t813) * t847 + t1073 * t813 * t1052 + t816 * t1006;
t894 = (t760 * t963 - t773 * t775) * t800;
t781 = t848 * t806 + t825 * t854;
t893 = (-t775 * t781 + t854 * t967) * t800;
t792 = t856 * qJ(2,1) - t825 * t850;
t765 = (t814 * t1051 - t792 * t817) * t849 + t817 * t954 + t814 * t1005;
t892 = (t761 * t960 - t765 * t776) * t803;
t771 = (t817 * t1051 + t792 * t814) * t849 + t1072 * t814 * t1050 + t817 * t1005;
t891 = (t761 * t961 - t771 * t776) * t803;
t782 = t850 * t807 + t825 * t856;
t890 = (-t776 * t782 + t856 * t966) * t803;
t869 = t867 / t866;
t737 = t887 - 0.2e1 * t1056;
t735 = t888 - 0.2e1 * t1057;
t733 = t889 - 0.2e1 * t1058;
t728 = (t863 + t870) * t746 - pkin(1) * t887;
t727 = (t862 + t870) * t745 - pkin(1) * t888;
t726 = (t861 + t870) * t744 - pkin(1) * t889;
t1 = [t936 + t938 + t940, 0, 0, (t737 * t1011 + t765 * t1044) * t801 + (t735 * t1012 + t767 * t1045) * t798 + (t733 * t1013 + t766 * t1046) * t795, t940 * t1060 + t938 * t1062 + t936 * t1064 + t826 * t898 + t830 * t895 + t834 * t892, t726 * t976 + t727 * t974 + t728 * t972 + (qJ(2,1) * t892 + t1018 * t765) * t834 + (qJ(2,2) * t895 + t1020 * t767) * t830 + (qJ(2,3) * t898 + t1022 * t766) * t826, t838 * t940 + t839 * t938 + t840 * t936 + (t851 * t904 + t853 * t902 + t855 * t900) * t867, -0.2e1 * t815 * t929 - 0.2e1 * t816 * t928 - 0.2e1 * t817 * t927 + (t900 * t1066 + t902 * t1067 + t904 * t1068) * t867, (-t812 * t992 - t813 * t991 - t814 * t990) * t867 + (-t1011 * t914 - t1012 * t915 - t1013 * t916) * t868, (t744 * t812 + t745 * t813 + t746 * t814) * t867, (t812 * t985 + t813 * t982 + t814 * t979) * t869, (t1000 * t1011 + t765 * t918) * t801 + (t1002 * t1012 + t767 * t919) * t798 + (t1004 * t1013 + t766 * t920) * t795 + (t812 * t913 + t813 * t912 + t814 * t911) * t867, (t1011 * t999 + t765 * t905) * t801 + (t1001 * t1012 + t767 * t906) * t798 + (t1003 * t1013 + t766 * t907) * t795 + (-t812 * t930 - t813 * t931 - t814 * t932) * t867, 0; -t937 - t939 - t941, 0, 0, (-t737 * t1014 + t771 * t1044) * t801 + (-t735 * t1015 + t773 * t1045) * t798 + (-t733 * t1016 + t772 * t1046) * t795, t941 * t1061 + t939 * t1063 + t937 * t1065 + t826 * t897 + t830 * t894 + t834 * t891, -t726 * t977 - t727 * t975 - t728 * t973 + (qJ(2,1) * t891 + t1018 * t771) * t834 + (qJ(2,2) * t894 + t1020 * t773) * t830 + (qJ(2,3) * t897 + t1022 * t772) * t826, -t838 * t941 - t839 * t939 - t840 * t937 + (t851 * t903 + t853 * t901 + t855 * t899) * t867, 0.2e1 * t812 * t929 + 0.2e1 * t813 * t928 + 0.2e1 * t814 * t927 + (t899 * t1066 + t901 * t1067 + t903 * t1068) * t867, (-t815 * t992 - t816 * t991 - t817 * t990) * t867 + (t973 * t1069 + t975 * t1070 + t977 * t1071) * t868, (t744 * t815 + t745 * t816 + t746 * t817) * t867, (t815 * t985 + t816 * t982 + t817 * t979) * t869, (-t1000 * t1014 + t771 * t918) * t801 + (-t1002 * t1015 + t773 * t919) * t798 + (-t1004 * t1016 + t772 * t920) * t795 + (t815 * t913 + t816 * t912 + t817 * t911) * t867, (-t1014 * t999 + t771 * t905) * t801 + (-t1001 * t1015 + t773 * t906) * t798 + (-t1003 * t1016 + t772 * t907) * t795 + (-t815 * t930 - t816 * t931 - t817 * t932) * t867, 0; t969 + t970 + t971, 0, 0, (t737 * t856 + t746 * t782) * t801 + (t735 * t854 + t745 * t781) * t798 + (t733 * t852 + t744 * t780) * t795, t971 * t1060 + t970 * t1062 + t969 * t1064 + t890 + t893 + t896, (t856 * t728 + t782 * t736) * t801 + (t854 * t727 + t781 * t734) * t798 + (t852 * t726 + t780 * t732) * t795 + qJ(2,3) * t896 + qJ(2,2) * t893 + qJ(2,1) * t890, t838 * t971 + t839 * t970 + t840 * t969 + 0.2e1 * (t851 * t947 + t853 * t946 + t855 * t945) * t867, -0.2e1 * t852 * t953 - 0.2e1 * t854 * t952 - 0.2e1 * t856 * t951 + 0.2e1 * (t945 * t1066 + t946 * t1067 + t947 * t1068) * t867, (-t852 * t916 - t854 * t915 - t856 * t914) * t868, 0, 0, (t1000 * t856 + (t770 * t849 - t942) * t782) * t801 + (t1002 * t854 + (t769 * t847 - t943) * t781) * t798 + (t1004 * t852 + (t768 * t845 - t944) * t780) * t795, (t999 * t856 + (t770 + t980) * t855 * t782) * t801 + (t1001 * t854 + (t769 + t983) * t853 * t781) * t798 + (t1003 * t852 + (t768 + t986) * t851 * t780) * t795, 0;];
tau_reg  = t1;
