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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 09:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR7V1G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 08:51:50
% EndTime: 2020-08-07 08:51:53
% DurationCPUTime: 2.64s
% Computational Cost: add. (6204->363), mult. (9981->543), div. (720->14), fcn. (5526->78), ass. (0->268)
t900 = sin(qJ(1,3));
t921 = -pkin(5) - pkin(4);
t857 = t900 * t921;
t909 = cos(qJ(1,3));
t907 = cos(qJ(3,3));
t1037 = t907 * pkin(2);
t852 = pkin(1) + t1037;
t908 = cos(qJ(2,3));
t898 = sin(qJ(3,3));
t899 = sin(qJ(2,3));
t988 = t898 * t899;
t945 = pkin(2) * t988 - t852 * t908;
t1052 = t909 * t945 + t857;
t903 = sin(qJ(1,2));
t858 = t903 * t921;
t912 = cos(qJ(1,2));
t910 = cos(qJ(3,2));
t1036 = t910 * pkin(2);
t854 = pkin(1) + t1036;
t911 = cos(qJ(2,2));
t901 = sin(qJ(3,2));
t902 = sin(qJ(2,2));
t987 = t901 * t902;
t944 = pkin(2) * t987 - t854 * t911;
t1051 = t912 * t944 + t858;
t906 = sin(qJ(1,1));
t859 = t906 * t921;
t915 = cos(qJ(1,1));
t913 = cos(qJ(3,1));
t1035 = t913 * pkin(2);
t856 = pkin(1) + t1035;
t914 = cos(qJ(2,1));
t904 = sin(qJ(3,1));
t905 = sin(qJ(2,1));
t986 = t904 * t905;
t943 = pkin(2) * t986 - t856 * t914;
t1050 = t915 * t943 + t859;
t1049 = m(2) / 0.2e1;
t1048 = m(3) / 0.2e1;
t1047 = Icges(2,2) / 0.2e1;
t1046 = Icges(3,2) / 0.2e1;
t1045 = 0.2e1 * pkin(2);
t1044 = 0.2e1 * t921;
t1043 = pkin(1) * m(3);
t1042 = m(2) * rSges(2,3);
t1041 = m(3) * rSges(3,2);
t932 = rSges(2,2) ^ 2;
t934 = rSges(2,1) ^ 2;
t936 = pkin(1) ^ 2;
t1040 = t936 * t1048 + (-t932 + t934) * t1049 + t1047 - Icges(2,1) / 0.2e1;
t931 = rSges(3,2) ^ 2;
t933 = rSges(3,1) ^ 2;
t1039 = (-t931 + t933) * t1048 - Icges(3,1) / 0.2e1 + t1046;
t916 = rSges(3,3) + pkin(4);
t1038 = m(3) * t916;
t846 = rSges(3,2) * t1038 - Icges(3,6);
t847 = rSges(3,1) * t1038 - Icges(3,5);
t939 = -rSges(2,1) * t1042 - pkin(1) * t1038 + Icges(2,5);
t979 = rSges(2,2) * t1042 - Icges(2,6);
t794 = (t846 * t898 - t847 * t907 + t939) * t899 - (t846 * t907 + t847 * t898 + t979) * t908;
t821 = 0.1e1 / t945;
t1034 = t794 * t821;
t795 = (t846 * t901 - t847 * t910 + t939) * t902 - (t846 * t910 + t847 * t901 + t979) * t911;
t822 = 0.1e1 / t944;
t1033 = t795 * t822;
t796 = (t846 * t904 - t847 * t913 + t939) * t905 - (t846 * t913 + t847 * t904 + t979) * t914;
t823 = 0.1e1 / t943;
t1032 = t796 * t823;
t824 = pkin(2) * t898 * t908 + t852 * t899;
t895 = legFrame(3,2);
t873 = sin(t895);
t876 = cos(t895);
t797 = -t1052 * t873 - t876 * t824;
t880 = 0.1e1 / t898;
t1031 = t797 * t880;
t798 = t1052 * t876 - t873 * t824;
t1030 = t798 * t880;
t825 = pkin(2) * t901 * t911 + t854 * t902;
t896 = legFrame(2,2);
t874 = sin(t896);
t877 = cos(t896);
t799 = -t1051 * t874 - t877 * t825;
t881 = 0.1e1 / t901;
t1029 = t799 * t881;
t800 = t1051 * t877 - t874 * t825;
t1028 = t800 * t881;
t826 = pkin(2) * t904 * t914 + t856 * t905;
t897 = legFrame(1,2);
t875 = sin(t897);
t878 = cos(t897);
t801 = -t1050 * t875 - t878 * t826;
t882 = 0.1e1 / t904;
t1027 = t801 * t882;
t802 = t1050 * t878 - t875 * t826;
t1026 = t802 * t882;
t976 = rSges(3,1) * t1043;
t848 = t907 * t976;
t978 = t931 + t933;
t953 = t936 + t978;
t977 = t932 + t934;
t946 = m(2) * t977 + m(3) * t953 + Icges(2,3) + Icges(3,3);
t975 = pkin(1) * t1041;
t949 = t898 * t975;
t803 = 0.2e1 * t848 + t946 - 0.2e1 * t949;
t1025 = t803 * t821;
t849 = t910 * t976;
t948 = t901 * t975;
t804 = 0.2e1 * t849 + t946 - 0.2e1 * t948;
t1024 = t804 * t822;
t850 = t913 * t976;
t947 = t904 * t975;
t805 = 0.2e1 * t850 + t946 - 0.2e1 * t947;
t1023 = t805 * t823;
t1022 = (-t900 * t945 + t909 * t921) * t880;
t1021 = (-t903 * t944 + t912 * t921) * t881;
t1020 = (-t906 * t943 + t915 * t921) * t882;
t890 = qJ(2,3) + qJ(3,3);
t867 = sin(t890);
t870 = cos(t890);
t809 = -t846 * t870 - t847 * t867;
t935 = 0.1e1 / pkin(2);
t1019 = t809 * t935;
t892 = qJ(2,2) + qJ(3,2);
t868 = sin(t892);
t871 = cos(t892);
t810 = -t846 * t871 - t847 * t868;
t1018 = t810 * t935;
t894 = qJ(2,1) + qJ(3,1);
t869 = sin(t894);
t872 = cos(t894);
t811 = -t846 * t872 - t847 * t869;
t1017 = t811 * t935;
t844 = m(3) * t978 + Icges(3,3);
t818 = t848 - t949 + t844;
t1016 = t818 * t821;
t1015 = t818 * t935;
t819 = t849 - t948 + t844;
t1014 = t819 * t822;
t1013 = t819 * t935;
t820 = t850 - t947 + t844;
t1012 = t820 * t823;
t1011 = t820 * t935;
t1010 = t821 * t880;
t1009 = t822 * t881;
t1008 = t823 * t882;
t834 = 0.1e1 / (pkin(1) * t908 + pkin(2) * t870);
t1007 = t834 * t900;
t1006 = t834 * t909;
t835 = 0.1e1 / (pkin(1) * t911 + pkin(2) * t871);
t1005 = t835 * t903;
t1004 = t835 * t912;
t836 = 0.1e1 / (pkin(1) * t914 + pkin(2) * t872);
t1003 = t836 * t906;
t1002 = t836 * t915;
t883 = t907 ^ 2;
t838 = pkin(1) * t907 + t1045 * t883 - pkin(2);
t1001 = t838 * t899;
t1000 = t838 * t909;
t885 = t910 ^ 2;
t839 = pkin(1) * t910 + t1045 * t885 - pkin(2);
t999 = t839 * t902;
t998 = t839 * t912;
t887 = t913 ^ 2;
t840 = pkin(1) * t913 + t1045 * t887 - pkin(2);
t997 = t840 * t905;
t996 = t840 * t915;
t995 = t844 * t935;
t994 = (pkin(1) + 0.2e1 * t1037) * t898;
t993 = (pkin(1) + 0.2e1 * t1036) * t901;
t992 = (pkin(1) + 0.2e1 * t1035) * t904;
t937 = 0.1e1 / pkin(1);
t991 = t880 * t937;
t990 = t881 * t937;
t989 = t882 * t937;
t985 = -qJ(3,1) + qJ(1,1);
t984 = qJ(3,1) + qJ(1,1);
t983 = -qJ(3,2) + qJ(1,2);
t982 = qJ(3,2) + qJ(1,2);
t981 = -qJ(3,3) + qJ(1,3);
t980 = qJ(3,3) + qJ(1,3);
t923 = 0.2e1 * qJ(2,3);
t889 = t923 + qJ(3,3);
t926 = 0.2e1 * qJ(2,2);
t891 = t926 + qJ(3,2);
t929 = 0.2e1 * qJ(2,1);
t893 = t929 + qJ(3,1);
t974 = t898 * t1037;
t973 = t901 * t1036;
t972 = t904 * t1035;
t815 = -t857 * t988 + (t883 - 0.1e1) * t909 * pkin(2);
t884 = t908 ^ 2;
t956 = t909 * t988;
t942 = pkin(1) * t956 + (t1045 * t956 + t857) * t907;
t782 = (-t1000 * t873 + t876 * t994) * t884 + (t1001 * t876 + t873 * t942) * t908 + t815 * t873 - t876 * t974;
t971 = t782 * t1010;
t783 = (t1000 * t876 + t873 * t994) * t884 + (t1001 * t873 - t876 * t942) * t908 - t815 * t876 - t873 * t974;
t970 = t783 * t1010;
t816 = -t858 * t987 + (t885 - 0.1e1) * t912 * pkin(2);
t886 = t911 ^ 2;
t955 = t912 * t987;
t941 = pkin(1) * t955 + (t1045 * t955 + t858) * t910;
t784 = (-t874 * t998 + t877 * t993) * t886 + (t874 * t941 + t877 * t999) * t911 + t816 * t874 - t877 * t973;
t969 = t784 * t1009;
t785 = (t874 * t993 + t877 * t998) * t886 + (t874 * t999 - t877 * t941) * t911 - t816 * t877 - t874 * t973;
t968 = t785 * t1009;
t817 = -t859 * t986 + (t887 - 0.1e1) * t915 * pkin(2);
t888 = t914 ^ 2;
t954 = t915 * t986;
t940 = pkin(1) * t954 + (t1045 * t954 + t859) * t913;
t786 = (-t875 * t996 + t878 * t992) * t888 + (t875 * t940 + t878 * t997) * t914 + t817 * t875 - t878 * t972;
t967 = t786 * t1008;
t787 = (t875 * t992 + t878 * t996) * t888 + (t875 * t997 - t878 * t940) * t914 - t817 * t878 - t875 * t972;
t966 = t787 * t1008;
t965 = t935 * t1022;
t964 = t935 * t1021;
t963 = t935 * t1020;
t962 = t873 * t1007;
t961 = t876 * t1007;
t960 = t874 * t1005;
t959 = t877 * t1005;
t958 = t875 * t1003;
t957 = t878 * t1003;
t922 = 0.2e1 * qJ(3,3);
t924 = -0.2e1 * qJ(2,3);
t952 = ((cos(qJ(2,3) - t981) + cos(qJ(2,3) + t980)) * t1044 + (sin(qJ(1,3) + t924 - 0.2e1 * qJ(3,3)) + sin(qJ(1,3) + t923 + t922) + 0.2e1 * t900) * pkin(2) + (sin(t924 + t981) + sin(qJ(1,3) + t889) + sin(t981) + sin(t980)) * pkin(1)) / ((-sin(t922 + qJ(2,3)) + t899) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - t867) * pkin(1)) / 0.2e1;
t925 = 0.2e1 * qJ(3,2);
t927 = -0.2e1 * qJ(2,2);
t951 = ((cos(qJ(2,2) - t983) + cos(qJ(2,2) + t982)) * t1044 + (sin(qJ(1,2) + t927 - 0.2e1 * qJ(3,2)) + sin(qJ(1,2) + t926 + t925) + 0.2e1 * t903) * pkin(2) + (sin(t927 + t983) + sin(qJ(1,2) + t891) + sin(t983) + sin(t982)) * pkin(1)) / ((-sin(t925 + qJ(2,2)) + t902) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - t868) * pkin(1)) / 0.2e1;
t928 = 0.2e1 * qJ(3,1);
t930 = -0.2e1 * qJ(2,1);
t950 = ((cos(qJ(2,1) - t985) + cos(qJ(2,1) + t984)) * t1044 + (sin(qJ(1,1) + t930 - 0.2e1 * qJ(3,1)) + sin(qJ(1,1) + t929 + t928) + 0.2e1 * t906) * pkin(2) + (sin(t930 + t985) + sin(qJ(1,1) + t893) + sin(t985) + sin(t984)) * pkin(1)) / ((-sin(t928 + qJ(2,1)) + t905) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - t869) * pkin(1)) / 0.2e1;
t938 = Icges(1,3) + (0.2e1 * t916 ^ 2 + t953) * t1048 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (0.2e1 * rSges(2,3) ^ 2 + t977) * t1049 + t1046 + t1047 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t866 = -m(2) * rSges(2,1) * rSges(2,2) + Icges(2,4);
t865 = -rSges(3,1) * t1041 + Icges(3,4);
t864 = 0.2e1 * t894;
t863 = 0.2e1 * t892;
t862 = 0.2e1 * t890;
t790 = cos(t864) * t1039 + cos(t929) * t1040 + t850 + t865 * sin(t864) + t866 * sin(t929) + (cos(t893) * rSges(3,1) + (-sin(t893) - t904) * rSges(3,2)) * t1043 + t938;
t789 = cos(t863) * t1039 + cos(t926) * t1040 + t849 + t865 * sin(t863) + t866 * sin(t926) + (cos(t891) * rSges(3,1) + (-sin(t891) - t901) * rSges(3,2)) * t1043 + t938;
t788 = cos(t862) * t1039 + cos(t923) * t1040 + t848 + t865 * sin(t862) + t866 * sin(t923) + (cos(t889) * rSges(3,1) + (-sin(t889) - t898) * rSges(3,2)) * t1043 + t938;
t781 = -t811 * t1002 + (t820 * t950 + t844 * t963) * t937;
t780 = -t810 * t1004 + (t819 * t951 + t844 * t964) * t937;
t779 = -t809 * t1006 + (t818 * t952 + t844 * t965) * t937;
t778 = -t796 * t1002 + (t805 * t950 + t820 * t963) * t937;
t777 = -t795 * t1004 + (t804 * t951 + t819 * t964) * t937;
t776 = -t794 * t1006 + (t803 * t952 + t818 * t965) * t937;
t775 = -t811 * t957 + (-t1012 * t787 + t802 * t995) * t989;
t774 = -t810 * t959 + (-t1014 * t785 + t800 * t995) * t990;
t773 = -t809 * t961 + (-t1016 * t783 + t798 * t995) * t991;
t772 = t811 * t958 + (-t1012 * t786 + t801 * t995) * t989;
t771 = t810 * t960 + (-t1014 * t784 + t799 * t995) * t990;
t770 = t809 * t962 + (-t1016 * t782 + t797 * t995) * t991;
t769 = -t796 * t957 + (t1011 * t802 - t1023 * t787) * t989;
t768 = -t795 * t959 + (t1013 * t800 - t1024 * t785) * t990;
t767 = -t794 * t961 + (t1015 * t798 - t1025 * t783) * t991;
t766 = t796 * t958 + (t1011 * t801 - t1023 * t786) * t989;
t765 = t795 * t960 + (t1013 * t799 - t1024 * t784) * t990;
t764 = t794 * t962 + (t1015 * t797 - t1025 * t782) * t991;
t763 = -t790 * t1002 + (t796 * t950 + t811 * t963) * t937;
t762 = -t789 * t1004 + (t795 * t951 + t810 * t964) * t937;
t761 = -t788 * t1006 + (t794 * t952 + t809 * t965) * t937;
t760 = -t790 * t957 + (t1017 * t802 - t1032 * t787) * t989;
t759 = -t789 * t959 + (t1018 * t800 - t1033 * t785) * t990;
t758 = -t788 * t961 + (t1019 * t798 - t1034 * t783) * t991;
t757 = t790 * t958 + (t1017 * t801 - t1032 * t786) * t989;
t756 = t789 * t960 + (t1018 * t799 - t1033 * t784) * t990;
t755 = t788 * t962 + (t1019 * t797 - t1034 * t782) * t991;
t1 = [-t758 * t961 - t759 * t959 - t760 * t957 + m(4) + (-t767 * t970 - t768 * t968 - t769 * t966 + (t1026 * t775 + t1028 * t774 + t1030 * t773) * t935) * t937, t758 * t962 + t759 * t960 + t760 * t958 + (-t767 * t971 - t768 * t969 - t769 * t967 + (t1027 * t775 + t1029 * t774 + t773 * t1031) * t935) * t937, -t758 * t1006 - t759 * t1004 - t760 * t1002 + (t769 * t950 + t768 * t951 + t767 * t952 + (t1020 * t775 + t1021 * t774 + t1022 * t773) * t935) * t937; -t755 * t961 - t756 * t959 - t757 * t957 + (-t764 * t970 - t765 * t968 - t766 * t966 + (t1026 * t772 + t1028 * t771 + t1030 * t770) * t935) * t937, t755 * t962 + t756 * t960 + t757 * t958 + m(4) + (-t764 * t971 - t765 * t969 - t766 * t967 + (t1027 * t772 + t1029 * t771 + t1031 * t770) * t935) * t937, -t755 * t1006 - t756 * t1004 - t757 * t1002 + (t766 * t950 + t765 * t951 + t764 * t952 + (t1020 * t772 + t1021 * t771 + t1022 * t770) * t935) * t937; -t761 * t961 - t762 * t959 - t763 * t957 + (-t776 * t970 - t777 * t968 - t778 * t966 + (t1026 * t781 + t1028 * t780 + t1030 * t779) * t935) * t937, t761 * t962 + t762 * t960 + t763 * t958 + (-t776 * t971 - t777 * t969 - t778 * t967 + (t1027 * t781 + t1029 * t780 + t1031 * t779) * t935) * t937, -t761 * t1006 - t762 * t1004 - t763 * t1002 + m(4) + (t778 * t950 + t777 * t951 + t776 * t952 + (t1020 * t781 + t1021 * t780 + t1022 * t779) * t935) * t937;];
MX  = t1;
