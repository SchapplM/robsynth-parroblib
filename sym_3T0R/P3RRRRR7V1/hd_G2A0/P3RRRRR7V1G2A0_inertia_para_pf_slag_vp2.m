% Calculate inertia matrix for parallel robot
% P3RRRRR7V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:52
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR7V1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:42:36
% EndTime: 2020-08-07 03:42:38
% DurationCPUTime: 2.56s
% Computational Cost: add. (5313->356), mult. (8649->541), div. (720->14), fcn. (5976->60), ass. (0->254)
t1053 = m(3) * pkin(4);
t1052 = 2 * pkin(2);
t1051 = -2 * Ifges(3,4);
t1050 = 4 * Ifges(3,4);
t921 = (-pkin(5) - pkin(4));
t1049 = 2 * t921;
t1048 = (pkin(1) * mrSges(3,1));
t1047 = (pkin(1) * mrSges(3,2));
t912 = cos(qJ(3,3));
t1046 = t912 * pkin(2);
t915 = cos(qJ(3,2));
t1045 = t915 * pkin(2);
t918 = cos(qJ(3,1));
t1044 = t918 * pkin(2);
t899 = Ifges(3,2) - Ifges(3,1);
t935 = 1 / pkin(2);
t1043 = Ifges(3,3) * t935;
t905 = sin(qJ(1,3));
t870 = pkin(1) + t1046;
t913 = cos(qJ(2,3));
t903 = sin(qJ(3,3));
t904 = sin(qJ(2,3));
t997 = t903 * t904;
t942 = pkin(2) * t997 - t870 * t913;
t914 = cos(qJ(1,3));
t988 = t914 * t921;
t824 = -t942 * t905 + t988;
t839 = t903 * pkin(2) * t913 + t904 * t870;
t900 = legFrame(3,2);
t877 = sin(t900);
t880 = cos(t900);
t812 = t824 * t877 - t839 * t880;
t887 = 0.1e1 / t903;
t1042 = t812 * t887;
t813 = -t824 * t880 - t839 * t877;
t1041 = t813 * t887;
t908 = sin(qJ(1,2));
t872 = pkin(1) + t1045;
t916 = cos(qJ(2,2));
t906 = sin(qJ(3,2));
t907 = sin(qJ(2,2));
t994 = t906 * t907;
t941 = pkin(2) * t994 - t872 * t916;
t917 = cos(qJ(1,2));
t987 = t917 * t921;
t825 = -t941 * t908 + t987;
t840 = t906 * pkin(2) * t916 + t907 * t872;
t901 = legFrame(2,2);
t878 = sin(t901);
t881 = cos(t901);
t814 = t825 * t878 - t840 * t881;
t888 = 0.1e1 / t906;
t1040 = t814 * t888;
t815 = -t825 * t881 - t840 * t878;
t1039 = t815 * t888;
t911 = sin(qJ(1,1));
t874 = pkin(1) + t1044;
t919 = cos(qJ(2,1));
t909 = sin(qJ(3,1));
t910 = sin(qJ(2,1));
t991 = t909 * t910;
t940 = pkin(2) * t991 - t874 * t919;
t920 = cos(qJ(1,1));
t986 = t920 * t921;
t826 = -t940 * t911 + t986;
t841 = t909 * pkin(2) * t919 + t910 * t874;
t902 = legFrame(1,2);
t879 = sin(t902);
t882 = cos(t902);
t816 = t826 * t879 - t841 * t882;
t889 = 0.1e1 / t909;
t1038 = t816 * t889;
t817 = -t826 * t882 - t841 * t879;
t1037 = t817 * t889;
t943 = Ifges(2,5) + (-mrSges(3,3) - t1053) * pkin(1);
t875 = pkin(4) * mrSges(3,2) - Ifges(3,6);
t876 = mrSges(3,1) * pkin(4) - Ifges(3,5);
t948 = -t875 * t912 - t903 * t876;
t979 = t903 * t875 - t876 * t912;
t818 = (t943 + t979) * t904 + (Ifges(2,6) + t948) * t913;
t836 = 0.1e1 / t942;
t1036 = t818 * t836;
t947 = -t875 * t915 - t906 * t876;
t978 = t906 * t875 - t876 * t915;
t819 = (t943 + t978) * t907 + (Ifges(2,6) + t947) * t916;
t837 = 0.1e1 / t941;
t1035 = t819 * t837;
t946 = -t875 * t918 - t909 * t876;
t977 = t909 * t875 - t876 * t918;
t820 = (t943 + t977) * t910 + (Ifges(2,6) + t946) * t919;
t838 = 0.1e1 / t940;
t1034 = t820 * t838;
t821 = t979 * t904 + t948 * t913;
t1033 = t821 * t935;
t822 = t978 * t907 + t947 * t916;
t1032 = t822 * t935;
t823 = t977 * t910 + t946 * t919;
t1031 = t823 * t935;
t1030 = (t905 * t921 + t942 * t914) * t887;
t1029 = (t908 * t921 + t941 * t917) * t888;
t1028 = (t911 * t921 + t940 * t920) * t889;
t975 = -2 * t1047;
t866 = t903 * t975;
t885 = 2 * t1048;
t925 = m(3) * pkin(1) ^ 2;
t970 = Ifges(2,3) + Ifges(3,3) + t925;
t845 = t885 * t912 + t866 + t970;
t1027 = t836 * t845;
t848 = Ifges(3,3) + (mrSges(3,1) * t912 - mrSges(3,2) * t903) * pkin(1);
t1026 = t836 * t848;
t1025 = t836 * t887;
t867 = t906 * t975;
t846 = t885 * t915 + t867 + t970;
t1024 = t837 * t846;
t849 = Ifges(3,3) + (mrSges(3,1) * t915 - mrSges(3,2) * t906) * pkin(1);
t1023 = t837 * t849;
t1022 = t837 * t888;
t868 = t909 * t975;
t847 = t885 * t918 + t868 + t970;
t1021 = t838 * t847;
t850 = Ifges(3,3) + (mrSges(3,1) * t918 - mrSges(3,2) * t909) * pkin(1);
t1020 = t838 * t850;
t1019 = t838 * t889;
t1018 = t848 * t935;
t1017 = t849 * t935;
t1016 = t850 * t935;
t896 = qJ(2,3) + qJ(3,3);
t851 = 0.1e1 / (pkin(2) * cos(t896) + t913 * pkin(1));
t1015 = t851 * t905;
t1014 = t851 * t914;
t897 = qJ(2,2) + qJ(3,2);
t852 = 0.1e1 / (pkin(2) * cos(t897) + t916 * pkin(1));
t1013 = t852 * t908;
t1012 = t852 * t917;
t898 = qJ(2,1) + qJ(3,1);
t853 = 0.1e1 / (pkin(2) * cos(t898) + t919 * pkin(1));
t1011 = t853 * t911;
t1010 = t853 * t920;
t890 = t912 ^ 2;
t854 = pkin(1) * t912 + t890 * t1052 - pkin(2);
t1009 = t854 * t905;
t892 = t915 ^ 2;
t855 = pkin(1) * t915 + t892 * t1052 - pkin(2);
t1008 = t855 * t908;
t894 = t918 ^ 2;
t856 = pkin(1) * t918 + t894 * t1052 - pkin(2);
t1007 = t856 * t911;
t1006 = (pkin(1) + 0.2e1 * t1046) * t903;
t1005 = (pkin(1) + 0.2e1 * t1045) * t906;
t1004 = (pkin(1) + 0.2e1 * t1044) * t909;
t936 = 0.1e1 / pkin(1);
t1003 = t887 * t936;
t1002 = t888 * t936;
t1001 = t889 * t936;
t1000 = t899 * t890;
t999 = t899 * t892;
t998 = t899 * t894;
t996 = t903 * t912;
t995 = t904 * t854;
t993 = t906 * t915;
t992 = t907 * t855;
t990 = t909 * t918;
t989 = t910 * t856;
t985 = -qJ(3,1) + qJ(1,1);
t984 = qJ(3,1) + qJ(1,1);
t983 = -qJ(3,2) + qJ(1,2);
t982 = qJ(3,2) + qJ(1,2);
t981 = -qJ(3,3) + qJ(1,3);
t980 = qJ(3,3) + qJ(1,3);
t976 = -Ifges(3,4) / 0.2e1 + Ifges(2,4) / 0.2e1;
t974 = pkin(2) * t996;
t973 = pkin(2) * t993;
t972 = pkin(2) * t990;
t971 = -t1048 / 0.2e1;
t833 = t988 * t997 + (t890 - 0.1e1) * t905 * pkin(2);
t891 = t913 ^ 2;
t954 = t905 * t997;
t939 = pkin(1) * t954 + (t954 * t1052 - t988) * t912;
t800 = (t880 * t1006 - t877 * t1009) * t891 + (t939 * t877 + t880 * t995) * t913 + t833 * t877 - t880 * t974;
t969 = t800 * t1025;
t801 = (t877 * t1006 + t880 * t1009) * t891 + (t877 * t995 - t939 * t880) * t913 - t833 * t880 - t877 * t974;
t968 = t801 * t1025;
t834 = t987 * t994 + (t892 - 0.1e1) * t908 * pkin(2);
t893 = t916 ^ 2;
t953 = t908 * t994;
t938 = pkin(1) * t953 + (t953 * t1052 - t987) * t915;
t802 = (t881 * t1005 - t878 * t1008) * t893 + (t938 * t878 + t881 * t992) * t916 + t834 * t878 - t881 * t973;
t967 = t802 * t1022;
t803 = (t878 * t1005 + t881 * t1008) * t893 + (t878 * t992 - t938 * t881) * t916 - t834 * t881 - t878 * t973;
t966 = t803 * t1022;
t835 = t986 * t991 + (t894 - 0.1e1) * t911 * pkin(2);
t895 = t919 ^ 2;
t952 = t911 * t991;
t937 = pkin(1) * t952 + (t952 * t1052 - t986) * t918;
t804 = (t882 * t1004 - t879 * t1007) * t895 + (t937 * t879 + t882 * t989) * t919 + t835 * t879 - t882 * t972;
t965 = t804 * t1019;
t805 = (t879 * t1004 + t882 * t1007) * t895 + (t879 * t989 - t937 * t882) * t919 - t835 * t882 - t879 * t972;
t964 = t805 * t1019;
t963 = t935 * t1030;
t962 = t935 * t1029;
t961 = t935 * t1028;
t960 = t877 * t1014;
t959 = t880 * t1014;
t958 = t878 * t1012;
t957 = t881 * t1012;
t956 = t879 * t1010;
t955 = t882 * t1010;
t926 = 0.2e1 * qJ(3,3);
t927 = 0.2e1 * qJ(2,3);
t928 = -0.2e1 * qJ(2,3);
t951 = ((-sin(qJ(2,3) - t981) + sin(qJ(2,3) + t980)) * t1049 + (-cos(t928 - 0.2e1 * qJ(3,3) + qJ(1,3)) - cos(t927 + t926 + qJ(1,3)) - 0.2e1 * t914) * pkin(2) + (-cos(t928 + t981) - cos(t927 + t980) - cos(t981) - cos(t980)) * pkin(1)) / ((-sin(t926 + qJ(2,3)) + t904) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - sin(t896)) * pkin(1)) / 0.2e1;
t929 = 0.2e1 * qJ(3,2);
t930 = 0.2e1 * qJ(2,2);
t931 = -0.2e1 * qJ(2,2);
t950 = ((-sin(qJ(2,2) - t983) + sin(qJ(2,2) + t982)) * t1049 + (-cos(t931 - 0.2e1 * qJ(3,2) + qJ(1,2)) - cos(t930 + t929 + qJ(1,2)) - 0.2e1 * t917) * pkin(2) + (-cos(t931 + t983) - cos(t930 + t982) - cos(t983) - cos(t982)) * pkin(1)) / ((-sin(t929 + qJ(2,2)) + t907) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - sin(t897)) * pkin(1)) / 0.2e1;
t932 = 0.2e1 * qJ(3,1);
t933 = 0.2e1 * qJ(2,1);
t934 = -0.2e1 * qJ(2,1);
t949 = ((-sin(qJ(2,1) - t985) + sin(qJ(2,1) + t984)) * t1049 + (-cos(t934 - 0.2e1 * qJ(3,1) + qJ(1,1)) - cos(t933 + t932 + qJ(1,1)) - 0.2e1 * t920) * pkin(2) + (-cos(t934 + t985) - cos(t933 + t984) - cos(t985) - cos(t984)) * pkin(1)) / ((-sin(t932 + qJ(2,1)) + t910) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - sin(t898)) * pkin(1)) / 0.2e1;
t945 = Ifges(2,2) + t925 - Ifges(2,1) - t899;
t944 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + (0.2e1 * mrSges(3,3) + t1053) * pkin(4);
t884 = -t1047 / 0.2e1;
t883 = -Ifges(3,2) / 0.2e1 + Ifges(3,1) / 0.2e1;
t808 = (0.2e1 * t998 + (t909 * t1050 + t885) * t918 + t868 + t945) * t895 + 0.4e1 * (Ifges(3,4) * t894 + (t883 * t909 + t884) * t918 + t909 * t971 + t976) * t910 * t919 - t998 + t990 * t1051 + t944;
t807 = (0.2e1 * t999 + (t906 * t1050 + t885) * t915 + t867 + t945) * t893 + 0.4e1 * (Ifges(3,4) * t892 + (t883 * t906 + t884) * t915 + t906 * t971 + t976) * t907 * t916 - t999 + t993 * t1051 + t944;
t806 = (0.2e1 * t1000 + (t903 * t1050 + t885) * t912 + t866 + t945) * t891 + 0.4e1 * (Ifges(3,4) * t890 + (t883 * t903 + t884) * t912 + t903 * t971 + t976) * t904 * t913 - t1000 + t996 * t1051 + t944;
t799 = -t823 * t1011 + (Ifges(3,3) * t961 + t850 * t949) * t936;
t798 = -t822 * t1013 + (Ifges(3,3) * t962 + t849 * t950) * t936;
t797 = -t821 * t1015 + (Ifges(3,3) * t963 + t848 * t951) * t936;
t796 = -t820 * t1011 + (t847 * t949 + t850 * t961) * t936;
t795 = -t819 * t1013 + (t846 * t950 + t849 * t962) * t936;
t794 = -t818 * t1015 + (t845 * t951 + t848 * t963) * t936;
t793 = t823 * t955 + (-t805 * t1020 + t817 * t1043) * t1001;
t792 = t822 * t957 + (-t803 * t1023 + t815 * t1043) * t1002;
t791 = t821 * t959 + (-t801 * t1026 + t813 * t1043) * t1003;
t790 = -t823 * t956 + (-t804 * t1020 + t816 * t1043) * t1001;
t789 = -t822 * t958 + (-t802 * t1023 + t814 * t1043) * t1002;
t788 = -t821 * t960 + (-t800 * t1026 + t812 * t1043) * t1003;
t787 = t820 * t955 + (t817 * t1016 - t805 * t1021) * t1001;
t786 = t819 * t957 + (t815 * t1017 - t803 * t1024) * t1002;
t785 = t818 * t959 + (t813 * t1018 - t801 * t1027) * t1003;
t784 = -t820 * t956 + (t816 * t1016 - t804 * t1021) * t1001;
t783 = -t819 * t958 + (t814 * t1017 - t802 * t1024) * t1002;
t782 = -t818 * t960 + (t812 * t1018 - t800 * t1027) * t1003;
t781 = -t808 * t1011 + (t820 * t949 + t823 * t961) * t936;
t780 = -t807 * t1013 + (t819 * t950 + t822 * t962) * t936;
t779 = -t806 * t1015 + (t818 * t951 + t821 * t963) * t936;
t778 = t808 * t955 + (t817 * t1031 - t805 * t1034) * t1001;
t777 = t807 * t957 + (t815 * t1032 - t803 * t1035) * t1002;
t776 = t806 * t959 + (t813 * t1033 - t801 * t1036) * t1003;
t775 = -t808 * t956 + (t816 * t1031 - t804 * t1034) * t1001;
t774 = -t807 * t958 + (t814 * t1032 - t802 * t1035) * t1002;
t773 = -t806 * t960 + (t812 * t1033 - t800 * t1036) * t1003;
t1 = [t776 * t959 + t777 * t957 + t778 * t955 + m(4) + (-t785 * t968 - t786 * t966 - t787 * t964 + (t793 * t1037 + t792 * t1039 + t791 * t1041) * t935) * t936, -t776 * t960 - t777 * t958 - t778 * t956 + (-t785 * t969 - t786 * t967 - t787 * t965 + (t793 * t1038 + t792 * t1040 + t791 * t1042) * t935) * t936, -t776 * t1015 - t777 * t1013 - t778 * t1011 + (t787 * t949 + t786 * t950 + t785 * t951 + (t1028 * t793 + t1029 * t792 + t1030 * t791) * t935) * t936; t773 * t959 + t774 * t957 + t775 * t955 + (-t782 * t968 - t783 * t966 - t784 * t964 + (t790 * t1037 + t789 * t1039 + t788 * t1041) * t935) * t936, -t773 * t960 - t774 * t958 - t775 * t956 + m(4) + (-t782 * t969 - t783 * t967 - t784 * t965 + (t1038 * t790 + t1040 * t789 + t1042 * t788) * t935) * t936, -t773 * t1015 - t774 * t1013 - t775 * t1011 + (t784 * t949 + t783 * t950 + t782 * t951 + (t1028 * t790 + t1029 * t789 + t1030 * t788) * t935) * t936; t779 * t959 + t780 * t957 + t781 * t955 + (-t794 * t968 - t795 * t966 - t796 * t964 + (t799 * t1037 + t798 * t1039 + t797 * t1041) * t935) * t936, -t779 * t960 - t780 * t958 - t781 * t956 + (-t794 * t969 - t795 * t967 - t796 * t965 + (t1038 * t799 + t1040 * t798 + t1042 * t797) * t935) * t936, -t779 * t1015 - t780 * t1013 - t781 * t1011 + m(4) + (t796 * t949 + t795 * t950 + t794 * t951 + (t1028 * t799 + t1029 * t798 + t1030 * t797) * t935) * t936;];
MX  = t1;
