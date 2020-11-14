% Calculate inertia matrix for parallel robot
% P3RRRRR7V1G3A0
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
% Datum: 2020-08-07 09:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR7V1G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 08:52:56
% EndTime: 2020-08-07 08:52:58
% DurationCPUTime: 2.50s
% Computational Cost: add. (5313->356), mult. (8577->541), div. (720->14), fcn. (5904->60), ass. (0->254)
t1059 = m(3) * pkin(4);
t911 = sin(qJ(1,3));
t927 = -pkin(5) - pkin(4);
t878 = t911 * t927;
t920 = cos(qJ(1,3));
t909 = sin(qJ(3,3));
t910 = sin(qJ(2,3));
t1000 = t909 * t910;
t918 = cos(qJ(3,3));
t1049 = t918 * pkin(2);
t873 = pkin(1) + t1049;
t919 = cos(qJ(2,3));
t948 = pkin(2) * t1000 - t873 * t919;
t1058 = t948 * t920 + t878;
t914 = sin(qJ(1,2));
t879 = t914 * t927;
t923 = cos(qJ(1,2));
t921 = cos(qJ(3,2));
t1048 = t921 * pkin(2);
t875 = pkin(1) + t1048;
t922 = cos(qJ(2,2));
t912 = sin(qJ(3,2));
t913 = sin(qJ(2,2));
t997 = t912 * t913;
t947 = pkin(2) * t997 - t875 * t922;
t1057 = t947 * t923 + t879;
t917 = sin(qJ(1,1));
t880 = t917 * t927;
t926 = cos(qJ(1,1));
t924 = cos(qJ(3,1));
t1047 = t924 * pkin(2);
t877 = pkin(1) + t1047;
t925 = cos(qJ(2,1));
t915 = sin(qJ(3,1));
t916 = sin(qJ(2,1));
t994 = t915 * t916;
t946 = pkin(2) * t994 - t877 * t925;
t1056 = t946 * t926 + t880;
t1055 = 0.2e1 * pkin(2);
t1054 = -2 * Ifges(3,4);
t1053 = 4 * Ifges(3,4);
t1052 = 0.2e1 * t927;
t1051 = pkin(1) * mrSges(3,1);
t1050 = pkin(1) * mrSges(3,2);
t905 = Ifges(3,2) - Ifges(3,1);
t941 = 0.1e1 / pkin(2);
t1046 = Ifges(3,3) * t941;
t839 = t909 * pkin(2) * t919 + t910 * t873;
t906 = legFrame(3,2);
t883 = sin(t906);
t886 = cos(t906);
t815 = -t1058 * t883 - t886 * t839;
t893 = 0.1e1 / t909;
t1045 = t815 * t893;
t816 = t1058 * t886 - t883 * t839;
t1044 = t816 * t893;
t840 = t912 * pkin(2) * t922 + t913 * t875;
t907 = legFrame(2,2);
t884 = sin(t907);
t887 = cos(t907);
t817 = -t1057 * t884 - t887 * t840;
t894 = 0.1e1 / t912;
t1043 = t817 * t894;
t818 = t1057 * t887 - t884 * t840;
t1042 = t818 * t894;
t841 = t915 * pkin(2) * t925 + t916 * t877;
t908 = legFrame(1,2);
t885 = sin(t908);
t888 = cos(t908);
t819 = -t1056 * t885 - t888 * t841;
t895 = 0.1e1 / t915;
t1041 = t819 * t895;
t820 = t1056 * t888 - t885 * t841;
t1040 = t820 * t895;
t949 = Ifges(2,5) + (-mrSges(3,3) - t1059) * pkin(1);
t881 = pkin(4) * mrSges(3,2) - Ifges(3,6);
t882 = mrSges(3,1) * pkin(4) - Ifges(3,5);
t954 = -t881 * t918 - t909 * t882;
t985 = t909 * t881 - t882 * t918;
t821 = (t949 + t985) * t910 + (Ifges(2,6) + t954) * t919;
t836 = 0.1e1 / t948;
t1039 = t821 * t836;
t953 = -t881 * t921 - t912 * t882;
t984 = t912 * t881 - t882 * t921;
t822 = (t949 + t984) * t913 + (Ifges(2,6) + t953) * t922;
t837 = 0.1e1 / t947;
t1038 = t822 * t837;
t952 = -t881 * t924 - t915 * t882;
t983 = t915 * t881 - t882 * t924;
t823 = (t949 + t983) * t916 + (Ifges(2,6) + t952) * t925;
t838 = 0.1e1 / t946;
t1037 = t823 * t838;
t824 = t985 * t910 + t954 * t919;
t1036 = t824 * t941;
t825 = t984 * t913 + t953 * t922;
t1035 = t825 * t941;
t826 = t983 * t916 + t952 * t925;
t1034 = t826 * t941;
t1033 = (-t911 * t948 + t920 * t927) * t893;
t1032 = (-t914 * t947 + t923 * t927) * t894;
t1031 = (-t917 * t946 + t926 * t927) * t895;
t981 = -0.2e1 * t1050;
t869 = t909 * t981;
t891 = 0.2e1 * t1051;
t931 = m(3) * pkin(1) ^ 2;
t976 = Ifges(2,3) + Ifges(3,3) + t931;
t848 = t918 * t891 + t869 + t976;
t1030 = t836 * t848;
t851 = Ifges(3,3) + (mrSges(3,1) * t918 - mrSges(3,2) * t909) * pkin(1);
t1029 = t836 * t851;
t1028 = t836 * t893;
t870 = t912 * t981;
t849 = t921 * t891 + t870 + t976;
t1027 = t837 * t849;
t852 = Ifges(3,3) + (mrSges(3,1) * t921 - mrSges(3,2) * t912) * pkin(1);
t1026 = t837 * t852;
t1025 = t837 * t894;
t871 = t915 * t981;
t850 = t924 * t891 + t871 + t976;
t1024 = t838 * t850;
t853 = Ifges(3,3) + (mrSges(3,1) * t924 - mrSges(3,2) * t915) * pkin(1);
t1023 = t838 * t853;
t1022 = t838 * t895;
t1021 = t851 * t941;
t1020 = t852 * t941;
t1019 = t853 * t941;
t902 = qJ(2,3) + qJ(3,3);
t854 = 0.1e1 / (pkin(2) * cos(t902) + t919 * pkin(1));
t1018 = t854 * t911;
t1017 = t854 * t920;
t903 = qJ(2,2) + qJ(3,2);
t855 = 0.1e1 / (pkin(2) * cos(t903) + t922 * pkin(1));
t1016 = t855 * t914;
t1015 = t855 * t923;
t904 = qJ(2,1) + qJ(3,1);
t856 = 0.1e1 / (pkin(2) * cos(t904) + t925 * pkin(1));
t1014 = t856 * t917;
t1013 = t856 * t926;
t896 = t918 ^ 2;
t857 = pkin(1) * t918 + t896 * t1055 - pkin(2);
t1012 = t857 * t920;
t898 = t921 ^ 2;
t858 = pkin(1) * t921 + t898 * t1055 - pkin(2);
t1011 = t858 * t923;
t900 = t924 ^ 2;
t859 = pkin(1) * t924 + t900 * t1055 - pkin(2);
t1010 = t859 * t926;
t1009 = (pkin(1) + 0.2e1 * t1049) * t909;
t1008 = (pkin(1) + 0.2e1 * t1048) * t912;
t1007 = (pkin(1) + 0.2e1 * t1047) * t915;
t942 = 0.1e1 / pkin(1);
t1006 = t893 * t942;
t1005 = t894 * t942;
t1004 = t895 * t942;
t1003 = t905 * t896;
t1002 = t905 * t898;
t1001 = t905 * t900;
t999 = t909 * t918;
t998 = t910 * t857;
t996 = t912 * t921;
t995 = t913 * t858;
t993 = t915 * t924;
t992 = t916 * t859;
t991 = -qJ(3,1) + qJ(1,1);
t990 = qJ(3,1) + qJ(1,1);
t989 = -qJ(3,2) + qJ(1,2);
t988 = qJ(3,2) + qJ(1,2);
t987 = -qJ(3,3) + qJ(1,3);
t986 = qJ(3,3) + qJ(1,3);
t982 = -Ifges(3,4) / 0.2e1 + Ifges(2,4) / 0.2e1;
t980 = pkin(2) * t999;
t979 = pkin(2) * t996;
t978 = pkin(2) * t993;
t977 = -t1051 / 0.2e1;
t833 = -t878 * t1000 + (t896 - 0.1e1) * t920 * pkin(2);
t897 = t919 ^ 2;
t960 = t920 * t1000;
t945 = pkin(1) * t960 + (t960 * t1055 + t878) * t918;
t803 = (t886 * t1009 - t883 * t1012) * t897 + (t945 * t883 + t886 * t998) * t919 + t833 * t883 - t886 * t980;
t975 = t803 * t1028;
t804 = (t883 * t1009 + t886 * t1012) * t897 + (t883 * t998 - t945 * t886) * t919 - t833 * t886 - t883 * t980;
t974 = t804 * t1028;
t834 = -t879 * t997 + (t898 - 0.1e1) * t923 * pkin(2);
t899 = t922 ^ 2;
t959 = t923 * t997;
t944 = pkin(1) * t959 + (t959 * t1055 + t879) * t921;
t805 = (t887 * t1008 - t884 * t1011) * t899 + (t944 * t884 + t887 * t995) * t922 + t834 * t884 - t887 * t979;
t973 = t805 * t1025;
t806 = (t884 * t1008 + t887 * t1011) * t899 + (t884 * t995 - t944 * t887) * t922 - t834 * t887 - t884 * t979;
t972 = t806 * t1025;
t835 = -t880 * t994 + (t900 - 0.1e1) * t926 * pkin(2);
t901 = t925 ^ 2;
t958 = t926 * t994;
t943 = pkin(1) * t958 + (t958 * t1055 + t880) * t924;
t807 = (t888 * t1007 - t885 * t1010) * t901 + (t943 * t885 + t888 * t992) * t925 + t835 * t885 - t888 * t978;
t971 = t807 * t1022;
t808 = (t885 * t1007 + t888 * t1010) * t901 + (t885 * t992 - t943 * t888) * t925 - t835 * t888 - t885 * t978;
t970 = t808 * t1022;
t969 = t941 * t1033;
t968 = t941 * t1032;
t967 = t941 * t1031;
t966 = t883 * t1018;
t965 = t886 * t1018;
t964 = t884 * t1016;
t963 = t887 * t1016;
t962 = t885 * t1014;
t961 = t888 * t1014;
t932 = 0.2e1 * qJ(3,3);
t933 = 0.2e1 * qJ(2,3);
t934 = -0.2e1 * qJ(2,3);
t957 = ((cos(qJ(2,3) - t987) + cos(qJ(2,3) + t986)) * t1052 + (sin(qJ(1,3) + t934 - 0.2e1 * qJ(3,3)) + sin(qJ(1,3) + t933 + t932) + 0.2e1 * t911) * pkin(2) + (sin(t934 + t987) + sin(t933 + t986) + sin(t987) + sin(t986)) * pkin(1)) / ((-sin(t932 + qJ(2,3)) + t910) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - sin(t902)) * pkin(1)) / 0.2e1;
t935 = 0.2e1 * qJ(3,2);
t936 = 0.2e1 * qJ(2,2);
t937 = -0.2e1 * qJ(2,2);
t956 = ((cos(qJ(2,2) - t989) + cos(qJ(2,2) + t988)) * t1052 + (sin(qJ(1,2) + t937 - 0.2e1 * qJ(3,2)) + sin(qJ(1,2) + t936 + t935) + 0.2e1 * t914) * pkin(2) + (sin(t937 + t989) + sin(t936 + t988) + sin(t989) + sin(t988)) * pkin(1)) / ((-sin(t935 + qJ(2,2)) + t913) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - sin(t903)) * pkin(1)) / 0.2e1;
t938 = 0.2e1 * qJ(3,1);
t939 = 0.2e1 * qJ(2,1);
t940 = -0.2e1 * qJ(2,1);
t955 = ((cos(qJ(2,1) - t991) + cos(qJ(2,1) + t990)) * t1052 + (sin(qJ(1,1) + t940 - 0.2e1 * qJ(3,1)) + sin(qJ(1,1) + t939 + t938) + 0.2e1 * t917) * pkin(2) + (sin(t940 + t991) + sin(t939 + t990) + sin(t991) + sin(t990)) * pkin(1)) / ((-sin(t938 + qJ(2,1)) + t916) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - sin(t904)) * pkin(1)) / 0.2e1;
t951 = Ifges(2,2) + t931 - Ifges(2,1) - t905;
t950 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + (0.2e1 * mrSges(3,3) + t1059) * pkin(4);
t890 = -t1050 / 0.2e1;
t889 = -Ifges(3,2) / 0.2e1 + Ifges(3,1) / 0.2e1;
t811 = (0.2e1 * t1001 + (t915 * t1053 + t891) * t924 + t871 + t951) * t901 + 0.4e1 * (Ifges(3,4) * t900 + (t889 * t915 + t890) * t924 + t915 * t977 + t982) * t916 * t925 - t1001 + t993 * t1054 + t950;
t810 = (0.2e1 * t1002 + (t912 * t1053 + t891) * t921 + t870 + t951) * t899 + 0.4e1 * (Ifges(3,4) * t898 + (t889 * t912 + t890) * t921 + t912 * t977 + t982) * t913 * t922 - t1002 + t996 * t1054 + t950;
t809 = (0.2e1 * t1003 + (t909 * t1053 + t891) * t918 + t869 + t951) * t897 + 0.4e1 * (Ifges(3,4) * t896 + (t889 * t909 + t890) * t918 + t909 * t977 + t982) * t910 * t919 - t1003 + t999 * t1054 + t950;
t802 = -t826 * t1013 + (Ifges(3,3) * t967 + t853 * t955) * t942;
t801 = -t825 * t1015 + (Ifges(3,3) * t968 + t852 * t956) * t942;
t800 = -t824 * t1017 + (Ifges(3,3) * t969 + t851 * t957) * t942;
t799 = -t823 * t1013 + (t850 * t955 + t853 * t967) * t942;
t798 = -t822 * t1015 + (t849 * t956 + t852 * t968) * t942;
t797 = -t821 * t1017 + (t848 * t957 + t851 * t969) * t942;
t796 = -t826 * t961 + (-t808 * t1023 + t820 * t1046) * t1004;
t795 = -t825 * t963 + (-t806 * t1026 + t818 * t1046) * t1005;
t794 = -t824 * t965 + (-t804 * t1029 + t816 * t1046) * t1006;
t793 = t826 * t962 + (-t807 * t1023 + t819 * t1046) * t1004;
t792 = t825 * t964 + (-t805 * t1026 + t817 * t1046) * t1005;
t791 = t824 * t966 + (-t803 * t1029 + t815 * t1046) * t1006;
t790 = -t823 * t961 + (t820 * t1019 - t808 * t1024) * t1004;
t789 = -t822 * t963 + (t818 * t1020 - t806 * t1027) * t1005;
t788 = -t821 * t965 + (t816 * t1021 - t804 * t1030) * t1006;
t787 = t823 * t962 + (t819 * t1019 - t807 * t1024) * t1004;
t786 = t822 * t964 + (t817 * t1020 - t805 * t1027) * t1005;
t785 = t821 * t966 + (t815 * t1021 - t803 * t1030) * t1006;
t784 = -t811 * t1013 + (t823 * t955 + t826 * t967) * t942;
t783 = -t810 * t1015 + (t822 * t956 + t825 * t968) * t942;
t782 = -t809 * t1017 + (t821 * t957 + t824 * t969) * t942;
t781 = -t811 * t961 + (t820 * t1034 - t808 * t1037) * t1004;
t780 = -t810 * t963 + (t818 * t1035 - t806 * t1038) * t1005;
t779 = -t809 * t965 + (t816 * t1036 - t804 * t1039) * t1006;
t778 = t811 * t962 + (t819 * t1034 - t807 * t1037) * t1004;
t777 = t810 * t964 + (t817 * t1035 - t805 * t1038) * t1005;
t776 = t809 * t966 + (t815 * t1036 - t803 * t1039) * t1006;
t1 = [-t779 * t965 - t780 * t963 - t781 * t961 + m(4) + (-t788 * t974 - t789 * t972 - t790 * t970 + (t796 * t1040 + t795 * t1042 + t794 * t1044) * t941) * t942, t779 * t966 + t780 * t964 + t781 * t962 + (-t788 * t975 - t789 * t973 - t790 * t971 + (t796 * t1041 + t795 * t1043 + t794 * t1045) * t941) * t942, -t779 * t1017 - t780 * t1015 - t781 * t1013 + (t790 * t955 + t789 * t956 + t788 * t957 + (t796 * t1031 + t795 * t1032 + t794 * t1033) * t941) * t942; -t776 * t965 - t777 * t963 - t778 * t961 + (-t785 * t974 - t786 * t972 - t787 * t970 + (t793 * t1040 + t792 * t1042 + t791 * t1044) * t941) * t942, t776 * t966 + t777 * t964 + t778 * t962 + m(4) + (-t785 * t975 - t786 * t973 - t787 * t971 + (t793 * t1041 + t792 * t1043 + t791 * t1045) * t941) * t942, -t776 * t1017 - t777 * t1015 - t778 * t1013 + (t787 * t955 + t786 * t956 + t785 * t957 + (t793 * t1031 + t792 * t1032 + t791 * t1033) * t941) * t942; -t782 * t965 - t783 * t963 - t784 * t961 + (-t797 * t974 - t798 * t972 - t799 * t970 + (t802 * t1040 + t801 * t1042 + t800 * t1044) * t941) * t942, t782 * t966 + t783 * t964 + t784 * t962 + (-t797 * t975 - t798 * t973 - t799 * t971 + (t802 * t1041 + t801 * t1043 + t800 * t1045) * t941) * t942, -t782 * t1017 - t783 * t1015 - t784 * t1013 + m(4) + (t799 * t955 + t798 * t956 + t797 * t957 + (t802 * t1031 + t801 * t1032 + t800 * t1033) * t941) * t942;];
MX  = t1;
