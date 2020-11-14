% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRRR12V1G3A0
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
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR12V1G3A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:28:23
% EndTime: 2020-08-06 18:28:27
% DurationCPUTime: 4.40s
% Computational Cost: add. (12977->334), mult. (21927->694), div. (2289->18), fcn. (15600->18), ass. (0->288)
t821 = sin(qJ(3,3));
t1032 = t821 * pkin(3);
t781 = qJ(2,3) + t1032;
t773 = 0.1e1 / t781 ^ 2;
t1039 = 0.2e1 * t773;
t1040 = 0.2e1 * qJ(2,3);
t818 = legFrame(3,2);
t788 = sin(t818);
t791 = cos(t818);
t834 = xDP(2);
t835 = xDP(1);
t766 = t788 * t835 + t791 * t834;
t827 = cos(qJ(3,3));
t1005 = t766 * t827;
t822 = sin(qJ(1,3));
t828 = cos(qJ(1,3));
t833 = xDP(3);
t899 = -t788 * t834 + t791 * t835;
t757 = -t822 * t833 + t899 * t828;
t802 = 0.1e1 / t821;
t1015 = t757 * t802;
t1036 = (-pkin(5) - pkin(6));
t801 = pkin(1) - t1036;
t1016 = t757 * t801;
t1017 = t757 * t773;
t814 = t827 ^ 2;
t1054 = -t814 + 0.1e1;
t794 = pkin(3) * t833;
t935 = pkin(3) * (t827 + 0.1e1) * (t827 - 0.1e1) * t822;
t986 = t827 * qJ(2,3);
t993 = t801 * t833;
t739 = ((qJ(2,3) * t833 + t899 * t801) * t828 + (t899 * qJ(2,3) - t993) * t822 + pkin(3) * t1005) * t821 - t899 * t935 + t766 * t986 + t1054 * t828 * t794;
t1023 = t739 * t802;
t772 = 0.1e1 / t781;
t774 = t772 * t773;
t724 = (t1005 * t1039 - t739 * t774) * t1015 - (-t1016 + t1023) * t772 * t1017;
t836 = pkin(1) + pkin(5);
t763 = t766 ^ 2;
t844 = 0.1e1 / pkin(3) ^ 2;
t1008 = t763 * t844;
t854 = t821 ^ 2;
t803 = 0.1e1 / t854;
t962 = t803 * t1008;
t914 = t836 * t962;
t968 = t739 * t1015;
t1057 = t968 * t1039 + t724 * t1040 + t914;
t823 = sin(qJ(3,2));
t1031 = t823 * pkin(3);
t782 = qJ(2,2) + t1031;
t776 = 0.1e1 / t782 ^ 2;
t1038 = 0.2e1 * t776;
t1042 = 0.2e1 * qJ(2,2);
t819 = legFrame(2,2);
t789 = sin(t819);
t792 = cos(t819);
t768 = t789 * t835 + t792 * t834;
t829 = cos(qJ(3,2));
t1004 = t768 * t829;
t824 = sin(qJ(1,2));
t830 = cos(qJ(1,2));
t898 = -t789 * t834 + t792 * t835;
t758 = -t824 * t833 + t898 * t830;
t806 = 0.1e1 / t823;
t1012 = t758 * t806;
t1013 = t758 * t801;
t1014 = t758 * t776;
t815 = t829 ^ 2;
t1053 = -t815 + 0.1e1;
t934 = pkin(3) * (t829 + 0.1e1) * (t829 - 0.1e1) * t824;
t985 = t829 * qJ(2,2);
t740 = ((qJ(2,2) * t833 + t898 * t801) * t830 + (t898 * qJ(2,2) - t993) * t824 + pkin(3) * t1004) * t823 - t898 * t934 + t768 * t985 + t1053 * t830 * t794;
t1022 = t740 * t806;
t775 = 0.1e1 / t782;
t777 = t775 * t776;
t725 = (t1004 * t1038 - t740 * t777) * t1012 - (-t1013 + t1022) * t775 * t1014;
t764 = t768 ^ 2;
t1007 = t764 * t844;
t857 = t823 ^ 2;
t807 = 0.1e1 / t857;
t959 = t807 * t1007;
t913 = t836 * t959;
t966 = t740 * t1012;
t1056 = t966 * t1038 + t725 * t1042 + t913;
t825 = sin(qJ(3,1));
t1030 = t825 * pkin(3);
t783 = qJ(2,1) + t1030;
t779 = 0.1e1 / t783 ^ 2;
t1037 = 0.2e1 * t779;
t1044 = 0.2e1 * qJ(2,1);
t820 = legFrame(1,2);
t790 = sin(t820);
t793 = cos(t820);
t770 = t790 * t835 + t793 * t834;
t831 = cos(qJ(3,1));
t1003 = t770 * t831;
t826 = sin(qJ(1,1));
t832 = cos(qJ(1,1));
t897 = -t790 * t834 + t793 * t835;
t759 = -t826 * t833 + t897 * t832;
t810 = 0.1e1 / t825;
t1009 = t759 * t810;
t1010 = t759 * t801;
t1011 = t759 * t779;
t816 = t831 ^ 2;
t1052 = -t816 + 0.1e1;
t933 = pkin(3) * (t831 + 0.1e1) * (t831 - 0.1e1) * t826;
t984 = t831 * qJ(2,1);
t741 = ((qJ(2,1) * t833 + t897 * t801) * t832 + (t897 * qJ(2,1) - t993) * t826 + pkin(3) * t1003) * t825 - t897 * t933 + t770 * t984 + t1052 * t832 * t794;
t1021 = t741 * t810;
t778 = 0.1e1 / t783;
t780 = t778 * t779;
t726 = (t1003 * t1037 - t741 * t780) * t1009 - (-t1010 + t1021) * t778 * t1011;
t765 = t770 ^ 2;
t1006 = t765 * t844;
t860 = t825 ^ 2;
t811 = 0.1e1 / t860;
t956 = t811 * t1006;
t912 = t836 * t956;
t964 = t741 * t1009;
t1055 = t964 * t1037 + t726 * t1044 + t912;
t804 = t802 * t803;
t1051 = t763 * (t804 * t814 + t802);
t808 = t806 * t807;
t1050 = t764 * (t808 * t815 + t806);
t812 = t810 * t811;
t1049 = t765 * (t812 * t816 + t810);
t1048 = (0.2e1 * t814 - 0.1e1) * t802;
t1047 = (0.2e1 * t815 - 0.1e1) * t806;
t1046 = (0.2e1 * t816 - 0.1e1) * t810;
t1041 = -0.2e1 * qJ(2,3);
t837 = qJ(2,3) ^ 2;
t842 = pkin(3) ^ 2;
t846 = pkin(1) ^ 2;
t893 = (2 * t1036 * pkin(1)) - (pkin(6) ^ 2) - t842 - t846 + ((-2 * pkin(6) - pkin(5)) * pkin(5));
t967 = t801 * t1023;
t843 = 0.1e1 / pkin(3);
t972 = qJ(2,3) * t766 * t843;
t992 = t802 * t827;
t865 = -((-t766 * t801 * t992 + (t967 + (t1032 * t1041 + t814 * t842 - t837 + t893) * t757) * t772) * t773 + t774 * t967) * t757 + ((t772 * t827 * t1016 + t766 * t802) * t821 + t802 * t972) * t803 * t772 * t766;
t1043 = -0.2e1 * qJ(2,2);
t838 = qJ(2,2) ^ 2;
t965 = t801 * t1022;
t973 = qJ(2,2) * t768 * t843;
t991 = t806 * t829;
t864 = -((-t768 * t801 * t991 + (t965 + (t1031 * t1043 + t815 * t842 - t838 + t893) * t758) * t775) * t776 + t777 * t965) * t758 + ((t775 * t829 * t1013 + t768 * t806) * t823 + t806 * t973) * t807 * t775 * t768;
t1045 = -0.2e1 * qJ(2,1);
t839 = qJ(2,1) ^ 2;
t963 = t801 * t1021;
t974 = qJ(2,1) * t770 * t843;
t990 = t810 * t831;
t863 = -((-t770 * t801 * t990 + (t963 + (t1030 * t1045 + t816 * t842 - t839 + t893) * t759) * t778) * t779 + t780 * t963) * t759 + ((t778 * t831 * t1010 + t770 * t810) * t825 + t810 * t974) * t811 * t778 * t770;
t1035 = pkin(1) * t724;
t1034 = pkin(1) * t725;
t1033 = pkin(1) * t726;
t1029 = t724 * t772;
t1028 = t724 * t802;
t1027 = t725 * t775;
t1026 = t725 * t806;
t1025 = t726 * t778;
t1024 = t726 * t810;
t754 = t757 ^ 2;
t1020 = t754 * t773;
t755 = t758 ^ 2;
t1019 = t755 * t776;
t756 = t759 ^ 2;
t1018 = t756 * t779;
t712 = t865 - t1035;
t1002 = t772 * t712;
t714 = t864 - t1034;
t1001 = t775 * t714;
t716 = t863 - t1033;
t1000 = t778 * t716;
t999 = t788 * t828;
t998 = t789 * t830;
t997 = t790 * t832;
t996 = t791 * t828;
t995 = t792 * t830;
t994 = t793 * t832;
t989 = t821 * t827;
t988 = t823 * t829;
t987 = t825 * t831;
t904 = t758 * t775 * t973;
t957 = t815 * t1007;
t922 = t808 * t957;
t983 = t1056 * t823 + t836 * t922 - 0.2e1 * t904 * t991;
t982 = 0.2e1 * t904 + (-t913 + t1056) * t829;
t903 = t759 * t778 * t974;
t954 = t816 * t1006;
t921 = t812 * t954;
t981 = t1055 * t825 + t836 * t921 - 0.2e1 * t903 * t990;
t980 = 0.2e1 * t903 + (-t912 + t1055) * t831;
t905 = t757 * t772 * t972;
t979 = 0.2e1 * t905 + (-t914 + t1057) * t827;
t960 = t814 * t1008;
t923 = t804 * t960;
t978 = t1057 * t821 + t836 * t923 - 0.2e1 * t905 * t992;
t971 = t724 * t992;
t970 = t725 * t991;
t969 = t726 * t990;
t805 = 0.1e1 / t854 ^ 2;
t961 = t763 * t805 * t827;
t809 = 0.1e1 / t857 ^ 2;
t958 = t764 * t809 * t829;
t813 = 0.1e1 / t860 ^ 2;
t955 = t765 * t813 * t831;
t953 = t772 * t999;
t952 = t772 * t996;
t951 = t775 * t998;
t950 = t775 * t995;
t949 = t778 * t997;
t948 = t778 * t994;
t947 = t822 * t1029;
t946 = t824 * t1027;
t945 = t826 * t1025;
t944 = -0.2e1 * t757 * t999;
t943 = 0.2e1 * t757 * t996;
t942 = -0.2e1 * t758 * t998;
t941 = 0.2e1 * t758 * t995;
t940 = -0.2e1 * t759 * t997;
t939 = 0.2e1 * t759 * t994;
t932 = t989 * t1029;
t931 = t988 * t1027;
t930 = t987 * t1025;
t926 = t766 * t822 * t1017;
t925 = t768 * t824 * t1014;
t924 = t770 * t826 * t1011;
t920 = t724 * t953;
t919 = t724 * t952;
t918 = t725 * t951;
t917 = t725 * t950;
t916 = t726 * t949;
t915 = t726 * t948;
t911 = qJ(2,1) * t1018 + t836 * t726 - t863;
t910 = qJ(2,2) * t1019 + t836 * t725 - t864;
t909 = qJ(2,3) * t1020 + t836 * t724 - t865;
t908 = t828 * t932;
t907 = t830 * t931;
t906 = t832 * t930;
t902 = t826 * qJ(2,1) + t801 * t832;
t901 = qJ(2,2) * t824 + t801 * t830;
t900 = qJ(2,3) * t822 + t801 * t828;
t751 = -t962 - t1020;
t896 = -t805 * t960 + t751;
t752 = -t959 - t1019;
t895 = -t809 * t957 + t752;
t753 = -t956 - t1018;
t894 = -t813 * t954 + t753;
t892 = t772 * t1051;
t891 = t775 * t1050;
t890 = t778 * t1049;
t889 = t909 * t992;
t888 = t910 * t991;
t887 = t911 * t990;
t883 = (t804 * t1008 + t751 * t802) * t827;
t882 = (t808 * t1007 + t752 * t806) * t829;
t881 = (t812 * t1006 + t753 * t810) * t831;
t880 = t773 * (-t754 * t788 + t766 * t943);
t879 = t773 * (-t754 * t791 + t766 * t944);
t878 = t776 * (-t755 * t789 + t768 * t941);
t877 = t776 * (-t755 * t792 + t768 * t942);
t876 = t779 * (-t756 * t790 + t770 * t939);
t875 = t779 * (-t756 * t793 + t770 * t940);
t745 = t900 * t791 * t821 + t788 * t986 + (t1054 * t791 * t822 + t788 * t989) * pkin(3);
t874 = (t739 * t943 - t745 * t754) * t774;
t748 = (t791 * pkin(3) * t827 - t900 * t788) * t821 + t788 * t935 + t791 * t986;
t873 = (t739 * t944 - t748 * t754) * t774;
t760 = t781 * t828 - t801 * t822;
t872 = (-t754 * t760 - 0.2e1 * t822 * t968) * t774;
t746 = t901 * t792 * t823 + t789 * t985 + (t1053 * t792 * t824 + t789 * t988) * pkin(3);
t871 = (t740 * t941 - t746 * t755) * t777;
t749 = (t792 * pkin(3) * t829 - t901 * t789) * t823 + t789 * t934 + t792 * t985;
t870 = (t740 * t942 - t749 * t755) * t777;
t761 = t782 * t830 - t801 * t824;
t869 = (-t755 * t761 - 0.2e1 * t824 * t966) * t777;
t747 = t902 * t793 * t825 + t790 * t984 + (t1052 * t793 * t826 + t790 * t987) * pkin(3);
t868 = (t741 * t939 - t747 * t756) * t780;
t750 = (t793 * pkin(3) * t831 - t902 * t790) * t825 + t790 * t933 + t793 * t984;
t867 = (t741 * t940 - t750 * t756) * t780;
t762 = t783 * t832 - t801 * t826;
t866 = (-t756 * t762 - 0.2e1 * t826 * t964) * t780;
t845 = t843 / t842;
t717 = t863 - 0.2e1 * t1033;
t715 = t864 - 0.2e1 * t1034;
t713 = t865 - 0.2e1 * t1035;
t708 = (t839 + t846) * t726 - pkin(1) * t863;
t707 = (t838 + t846) * t725 - pkin(1) * t864;
t706 = (t837 + t846) * t724 - pkin(1) * t865;
t1 = [t915 + t917 + t919, 0, 0, (t747 * t1024 + t717 * t994) * t778 + (t746 * t1026 + t715 * t995) * t775 + (t745 * t1028 + t713 * t996) * t772, t919 * t1040 + t917 * t1042 + t915 * t1044 + t802 * t874 + t806 * t871 + t810 * t868, t706 * t952 + t707 * t950 + t708 * t948 + (qJ(2,1) * t868 + t1000 * t747) * t810 + (qJ(2,2) * t871 + t1001 * t746) * t806 + (qJ(2,3) * t874 + t1002 * t745) * t802, t814 * t919 + t815 * t917 + t816 * t915 + (t827 * t880 + t829 * t878 + t831 * t876) * t843, -0.2e1 * t791 * t908 - 0.2e1 * t792 * t907 - 0.2e1 * t793 * t906 + (t876 * t1046 + t878 * t1047 + t880 * t1048) * t843, (-t788 * t971 - t789 * t970 - t790 * t969) * t843 + (-t948 * t1049 - t950 * t1050 - t952 * t1051) * t844, (t724 * t788 + t725 * t789 + t726 * t790) * t843, (t788 * t961 + t789 * t958 + t790 * t955) * t845, (t747 * t894 + t981 * t994) * t778 + (t746 * t895 + t983 * t995) * t775 + (t745 * t896 + t978 * t996) * t772 + (t788 * t889 + t789 * t888 + t790 * t887) * t843, (t747 * t881 + t980 * t994) * t778 + (t746 * t882 + t982 * t995) * t775 + (t745 * t883 + t979 * t996) * t772 + (-t788 * t909 - t789 * t910 - t790 * t911) * t843, 0; -t916 - t918 - t920, 0, 0, (t750 * t1024 - t717 * t997) * t778 + (t749 * t1026 - t715 * t998) * t775 + (t748 * t1028 - t713 * t999) * t772, t920 * t1041 + t918 * t1043 + t916 * t1045 + t802 * t873 + t806 * t870 + t810 * t867, -t706 * t953 - t707 * t951 - t708 * t949 + (qJ(2,1) * t867 + t1000 * t750) * t810 + (qJ(2,2) * t870 + t1001 * t749) * t806 + (qJ(2,3) * t873 + t1002 * t748) * t802, -t814 * t920 - t815 * t918 - t816 * t916 + (t827 * t879 + t829 * t877 + t831 * t875) * t843, 0.2e1 * t788 * t908 + 0.2e1 * t789 * t907 + 0.2e1 * t790 * t906 + (t875 * t1046 + t877 * t1047 + t879 * t1048) * t843, (-t791 * t971 - t792 * t970 - t793 * t969) * t843 + (t890 * t997 + t891 * t998 + t892 * t999) * t844, (t724 * t791 + t725 * t792 + t726 * t793) * t843, (t791 * t961 + t792 * t958 + t793 * t955) * t845, (t750 * t894 - t981 * t997) * t778 + (t749 * t895 - t983 * t998) * t775 + (t748 * t896 - t978 * t999) * t772 + (t791 * t889 + t792 * t888 + t793 * t887) * t843, (t750 * t881 - t980 * t997) * t778 + (t749 * t882 - t982 * t998) * t775 + (t748 * t883 - t979 * t999) * t772 + (-t791 * t909 - t792 * t910 - t793 * t911) * t843, 0; -t945 - t946 - t947, 0, 0, (-t717 * t826 + t726 * t762) * t778 + (-t715 * t824 + t725 * t761) * t775 + (-t713 * t822 + t724 * t760) * t772, t947 * t1041 + t946 * t1043 + t945 * t1045 + t866 + t869 + t872, (-t826 * t708 + t762 * t716) * t778 + (-t824 * t707 + t761 * t714) * t775 + (-t822 * t706 + t760 * t712) * t772 + qJ(2,3) * t872 + qJ(2,2) * t869 + qJ(2,1) * t866, -t814 * t947 - t815 * t946 - t816 * t945 + 0.2e1 * (-t827 * t926 - t829 * t925 - t831 * t924) * t843, 0.2e1 * t822 * t932 + 0.2e1 * t824 * t931 + 0.2e1 * t826 * t930 + 0.2e1 * (-t924 * t1046 - t925 * t1047 - t926 * t1048) * t843, (t822 * t892 + t824 * t891 + t826 * t890) * t844, 0, 0, (-t981 * t826 + (t753 * t825 - t921) * t762) * t778 + (-t983 * t824 + (t752 * t823 - t922) * t761) * t775 + (-t978 * t822 + (t751 * t821 - t923) * t760) * t772, (-t980 * t826 + (t753 + t956) * t831 * t762) * t778 + (-t982 * t824 + (t752 + t959) * t829 * t761) * t775 + (-t979 * t822 + (t751 + t962) * t827 * t760) * t772, 0;];
tau_reg  = t1;
