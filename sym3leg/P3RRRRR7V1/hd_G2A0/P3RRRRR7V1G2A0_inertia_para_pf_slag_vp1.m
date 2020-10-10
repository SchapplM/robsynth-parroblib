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
% Datum: 2020-08-07 03:52
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR7V1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:41:30
% EndTime: 2020-08-07 03:41:32
% DurationCPUTime: 2.35s
% Computational Cost: add. (6204->363), mult. (10053->543), div. (720->14), fcn. (5598->78), ass. (0->268)
t1046 = m(2) / 0.2e1;
t1045 = m(3) / 0.2e1;
t1044 = Icges(2,2) / 0.2e1;
t1043 = Icges(3,2) / 0.2e1;
t1042 = 2 * pkin(2);
t915 = (-pkin(5) - pkin(4));
t1041 = 2 * t915;
t1040 = pkin(1) * m(3);
t1039 = m(2) * rSges(2,3);
t1038 = m(3) * rSges(3,2);
t926 = rSges(2,2) ^ 2;
t928 = rSges(2,1) ^ 2;
t930 = pkin(1) ^ 2;
t1037 = t930 * t1045 + (-t926 + t928) * t1046 + t1044 - Icges(2,1) / 0.2e1;
t925 = rSges(3,2) ^ 2;
t927 = rSges(3,1) ^ 2;
t1036 = (-t925 + t927) * t1045 - Icges(3,1) / 0.2e1 + t1043;
t910 = rSges(3,3) + pkin(4);
t1035 = m(3) * t910;
t901 = cos(qJ(3,3));
t1034 = t901 * pkin(2);
t904 = cos(qJ(3,2));
t1033 = t904 * pkin(2);
t907 = cos(qJ(3,1));
t1032 = t907 * pkin(2);
t843 = rSges(3,2) * t1035 - Icges(3,6);
t844 = rSges(3,1) * t1035 - Icges(3,5);
t892 = sin(qJ(3,3));
t893 = sin(qJ(2,3));
t902 = cos(qJ(2,3));
t933 = -rSges(2,1) * t1039 - pkin(1) * t1035 + Icges(2,5);
t973 = rSges(2,2) * t1039 - Icges(2,6);
t791 = (t843 * t892 - t844 * t901 + t933) * t893 - (t843 * t901 + t844 * t892 + t973) * t902;
t849 = pkin(1) + t1034;
t985 = t892 * t893;
t939 = pkin(2) * t985 - t849 * t902;
t821 = 0.1e1 / t939;
t1031 = t791 * t821;
t895 = sin(qJ(3,2));
t896 = sin(qJ(2,2));
t905 = cos(qJ(2,2));
t792 = (t843 * t895 - t844 * t904 + t933) * t896 - (t843 * t904 + t844 * t895 + t973) * t905;
t851 = pkin(1) + t1033;
t984 = t895 * t896;
t938 = pkin(2) * t984 - t851 * t905;
t822 = 0.1e1 / t938;
t1030 = t792 * t822;
t898 = sin(qJ(3,1));
t899 = sin(qJ(2,1));
t908 = cos(qJ(2,1));
t793 = (t843 * t898 - t844 * t907 + t933) * t899 - (t843 * t907 + t844 * t898 + t973) * t908;
t853 = pkin(1) + t1032;
t983 = t898 * t899;
t937 = pkin(2) * t983 - t853 * t908;
t823 = 0.1e1 / t937;
t1029 = t793 * t823;
t894 = sin(qJ(1,3));
t903 = cos(qJ(1,3));
t982 = t903 * t915;
t803 = -t939 * t894 + t982;
t824 = t892 * pkin(2) * t902 + t893 * t849;
t889 = legFrame(3,2);
t867 = sin(t889);
t870 = cos(t889);
t794 = t803 * t867 - t824 * t870;
t874 = 0.1e1 / t892;
t1028 = t794 * t874;
t795 = -t803 * t870 - t824 * t867;
t1027 = t795 * t874;
t897 = sin(qJ(1,2));
t906 = cos(qJ(1,2));
t981 = t906 * t915;
t804 = -t938 * t897 + t981;
t825 = t895 * pkin(2) * t905 + t896 * t851;
t890 = legFrame(2,2);
t868 = sin(t890);
t871 = cos(t890);
t796 = t804 * t868 - t825 * t871;
t875 = 0.1e1 / t895;
t1026 = t796 * t875;
t797 = -t804 * t871 - t825 * t868;
t1025 = t797 * t875;
t900 = sin(qJ(1,1));
t909 = cos(qJ(1,1));
t980 = t909 * t915;
t805 = -t937 * t900 + t980;
t826 = t898 * pkin(2) * t908 + t899 * t853;
t891 = legFrame(1,2);
t869 = sin(t891);
t872 = cos(t891);
t798 = t805 * t869 - t826 * t872;
t876 = 0.1e1 / t898;
t1024 = t798 * t876;
t799 = -t805 * t872 - t826 * t869;
t1023 = t799 * t876;
t970 = rSges(3,1) * t1040;
t845 = t901 * t970;
t972 = t925 + t927;
t947 = t930 + t972;
t971 = t926 + t928;
t940 = t971 * m(2) + t947 * m(3) + Icges(2,3) + Icges(3,3);
t969 = pkin(1) * t1038;
t943 = t892 * t969;
t800 = 0.2e1 * t845 + t940 - 0.2e1 * t943;
t1022 = t800 * t821;
t846 = t904 * t970;
t942 = t895 * t969;
t801 = 0.2e1 * t846 + t940 - 0.2e1 * t942;
t1021 = t801 * t822;
t847 = t907 * t970;
t941 = t898 * t969;
t802 = 0.2e1 * t847 + t940 - 0.2e1 * t941;
t1020 = t802 * t823;
t1019 = (t894 * t915 + t939 * t903) * t874;
t1018 = (t897 * t915 + t938 * t906) * t875;
t1017 = (t900 * t915 + t937 * t909) * t876;
t884 = qJ(2,3) + qJ(3,3);
t861 = sin(t884);
t864 = cos(t884);
t809 = -t843 * t864 - t844 * t861;
t929 = 1 / pkin(2);
t1016 = t809 * t929;
t886 = qJ(2,2) + qJ(3,2);
t862 = sin(t886);
t865 = cos(t886);
t810 = -t843 * t865 - t844 * t862;
t1015 = t810 * t929;
t888 = qJ(2,1) + qJ(3,1);
t863 = sin(t888);
t866 = cos(t888);
t811 = -t843 * t866 - t844 * t863;
t1014 = t811 * t929;
t841 = t972 * m(3) + Icges(3,3);
t818 = t845 - t943 + t841;
t1013 = t818 * t821;
t1012 = t818 * t929;
t819 = t846 - t942 + t841;
t1011 = t819 * t822;
t1010 = t819 * t929;
t820 = t847 - t941 + t841;
t1009 = t820 * t823;
t1008 = t820 * t929;
t1007 = t821 * t874;
t1006 = t822 * t875;
t1005 = t823 * t876;
t831 = 0.1e1 / (t902 * pkin(1) + pkin(2) * t864);
t1004 = t831 * t894;
t1003 = t831 * t903;
t832 = 0.1e1 / (t905 * pkin(1) + pkin(2) * t865);
t1002 = t832 * t897;
t1001 = t832 * t906;
t833 = 0.1e1 / (t908 * pkin(1) + pkin(2) * t866);
t1000 = t833 * t900;
t999 = t833 * t909;
t877 = t901 ^ 2;
t835 = pkin(1) * t901 + t877 * t1042 - pkin(2);
t998 = t835 * t893;
t997 = t835 * t894;
t879 = t904 ^ 2;
t836 = pkin(1) * t904 + t879 * t1042 - pkin(2);
t996 = t836 * t896;
t995 = t836 * t897;
t881 = t907 ^ 2;
t837 = pkin(1) * t907 + t881 * t1042 - pkin(2);
t994 = t837 * t899;
t993 = t837 * t900;
t992 = t841 * t929;
t991 = (pkin(1) + 0.2e1 * t1034) * t892;
t990 = (pkin(1) + 0.2e1 * t1033) * t895;
t989 = (pkin(1) + 0.2e1 * t1032) * t898;
t931 = 0.1e1 / pkin(1);
t988 = t874 * t931;
t987 = t875 * t931;
t986 = t876 * t931;
t979 = -qJ(3,1) + qJ(1,1);
t978 = qJ(3,1) + qJ(1,1);
t977 = -qJ(3,2) + qJ(1,2);
t976 = qJ(3,2) + qJ(1,2);
t975 = -qJ(3,3) + qJ(1,3);
t974 = qJ(3,3) + qJ(1,3);
t917 = 0.2e1 * qJ(2,3);
t883 = t917 + qJ(3,3);
t920 = 0.2e1 * qJ(2,2);
t885 = t920 + qJ(3,2);
t923 = 0.2e1 * qJ(2,1);
t887 = t923 + qJ(3,1);
t968 = t892 * t1034;
t967 = t895 * t1033;
t966 = t898 * t1032;
t815 = t982 * t985 + (t877 - 0.1e1) * t894 * pkin(2);
t878 = t902 ^ 2;
t950 = t894 * t985;
t936 = pkin(1) * t950 + (t950 * t1042 - t982) * t901;
t779 = (-t867 * t997 + t870 * t991) * t878 + (t936 * t867 + t870 * t998) * t902 + t815 * t867 - t870 * t968;
t965 = t779 * t1007;
t780 = (t867 * t991 + t870 * t997) * t878 + (t867 * t998 - t936 * t870) * t902 - t815 * t870 - t867 * t968;
t964 = t780 * t1007;
t816 = t981 * t984 + (t879 - 0.1e1) * t897 * pkin(2);
t880 = t905 ^ 2;
t949 = t897 * t984;
t935 = pkin(1) * t949 + (t949 * t1042 - t981) * t904;
t781 = (-t868 * t995 + t871 * t990) * t880 + (t935 * t868 + t871 * t996) * t905 + t816 * t868 - t871 * t967;
t963 = t781 * t1006;
t782 = (t868 * t990 + t871 * t995) * t880 + (t868 * t996 - t935 * t871) * t905 - t816 * t871 - t868 * t967;
t962 = t782 * t1006;
t817 = t980 * t983 + (t881 - 0.1e1) * t900 * pkin(2);
t882 = t908 ^ 2;
t948 = t900 * t983;
t934 = pkin(1) * t948 + (t948 * t1042 - t980) * t907;
t783 = (-t869 * t993 + t872 * t989) * t882 + (t934 * t869 + t872 * t994) * t908 + t817 * t869 - t872 * t966;
t961 = t783 * t1005;
t784 = (t869 * t989 + t872 * t993) * t882 + (t869 * t994 - t934 * t872) * t908 - t817 * t872 - t869 * t966;
t960 = t784 * t1005;
t959 = t929 * t1019;
t958 = t929 * t1018;
t957 = t929 * t1017;
t956 = t867 * t1003;
t955 = t870 * t1003;
t954 = t868 * t1001;
t953 = t871 * t1001;
t952 = t869 * t999;
t951 = t872 * t999;
t916 = 0.2e1 * qJ(3,3);
t918 = -0.2e1 * qJ(2,3);
t946 = ((-sin(qJ(2,3) - t975) + sin(qJ(2,3) + t974)) * t1041 + (-cos(t918 - 0.2e1 * qJ(3,3) + qJ(1,3)) - cos(t917 + t916 + qJ(1,3)) - 0.2e1 * t903) * pkin(2) + (-cos(t918 + t975) - cos(qJ(1,3) + t883) - cos(t975) - cos(t974)) * pkin(1)) / ((-sin(t916 + qJ(2,3)) + t893) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - t861) * pkin(1)) / 0.2e1;
t919 = 0.2e1 * qJ(3,2);
t921 = -0.2e1 * qJ(2,2);
t945 = ((-sin(qJ(2,2) - t977) + sin(qJ(2,2) + t976)) * t1041 + (-cos(t921 - 0.2e1 * qJ(3,2) + qJ(1,2)) - cos(t920 + t919 + qJ(1,2)) - 0.2e1 * t906) * pkin(2) + (-cos(t921 + t977) - cos(qJ(1,2) + t885) - cos(t977) - cos(t976)) * pkin(1)) / ((-sin(t919 + qJ(2,2)) + t896) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - t862) * pkin(1)) / 0.2e1;
t922 = 0.2e1 * qJ(3,1);
t924 = -0.2e1 * qJ(2,1);
t944 = ((-sin(qJ(2,1) - t979) + sin(qJ(2,1) + t978)) * t1041 + (-cos(t924 - 0.2e1 * qJ(3,1) + qJ(1,1)) - cos(t923 + t922 + qJ(1,1)) - 0.2e1 * t909) * pkin(2) + (-cos(t924 + t979) - cos(qJ(1,1) + t887) - cos(t979) - cos(t978)) * pkin(1)) / ((-sin(t922 + qJ(2,1)) + t899) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - t863) * pkin(1)) / 0.2e1;
t932 = Icges(1,3) + (0.2e1 * t910 ^ 2 + t947) * t1045 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (0.2e1 * rSges(2,3) ^ 2 + t971) * t1046 + t1043 + t1044 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t860 = -m(2) * rSges(2,1) * rSges(2,2) + Icges(2,4);
t859 = -rSges(3,1) * t1038 + Icges(3,4);
t858 = 0.2e1 * t888;
t857 = 0.2e1 * t886;
t856 = 0.2e1 * t884;
t787 = cos(t858) * t1036 + cos(t923) * t1037 + t847 + t859 * sin(t858) + t860 * sin(t923) + (cos(t887) * rSges(3,1) + (-sin(t887) - t898) * rSges(3,2)) * t1040 + t932;
t786 = cos(t857) * t1036 + cos(t920) * t1037 + t846 + t859 * sin(t857) + t860 * sin(t920) + (cos(t885) * rSges(3,1) + (-sin(t885) - t895) * rSges(3,2)) * t1040 + t932;
t785 = cos(t856) * t1036 + cos(t917) * t1037 + t845 + t859 * sin(t856) + t860 * sin(t917) + (cos(t883) * rSges(3,1) + (-sin(t883) - t892) * rSges(3,2)) * t1040 + t932;
t778 = -t811 * t1000 + (t820 * t944 + t841 * t957) * t931;
t777 = -t810 * t1002 + (t819 * t945 + t841 * t958) * t931;
t776 = -t809 * t1004 + (t818 * t946 + t841 * t959) * t931;
t775 = -t793 * t1000 + (t802 * t944 + t820 * t957) * t931;
t774 = -t792 * t1002 + (t801 * t945 + t819 * t958) * t931;
t773 = -t791 * t1004 + (t800 * t946 + t818 * t959) * t931;
t772 = t811 * t951 + (-t784 * t1009 + t799 * t992) * t986;
t771 = t810 * t953 + (-t782 * t1011 + t797 * t992) * t987;
t770 = t809 * t955 + (-t780 * t1013 + t795 * t992) * t988;
t769 = -t811 * t952 + (-t783 * t1009 + t798 * t992) * t986;
t768 = -t810 * t954 + (-t781 * t1011 + t796 * t992) * t987;
t767 = -t809 * t956 + (-t779 * t1013 + t794 * t992) * t988;
t766 = t793 * t951 + (t799 * t1008 - t784 * t1020) * t986;
t765 = t792 * t953 + (t797 * t1010 - t782 * t1021) * t987;
t764 = t791 * t955 + (t795 * t1012 - t780 * t1022) * t988;
t763 = -t793 * t952 + (t798 * t1008 - t783 * t1020) * t986;
t762 = -t792 * t954 + (t796 * t1010 - t781 * t1021) * t987;
t761 = -t791 * t956 + (t794 * t1012 - t779 * t1022) * t988;
t760 = -t787 * t1000 + (t793 * t944 + t811 * t957) * t931;
t759 = -t786 * t1002 + (t792 * t945 + t810 * t958) * t931;
t758 = -t785 * t1004 + (t791 * t946 + t809 * t959) * t931;
t757 = t787 * t951 + (t799 * t1014 - t784 * t1029) * t986;
t756 = t786 * t953 + (t797 * t1015 - t782 * t1030) * t987;
t755 = t785 * t955 + (t795 * t1016 - t780 * t1031) * t988;
t754 = -t787 * t952 + (t798 * t1014 - t783 * t1029) * t986;
t753 = -t786 * t954 + (t796 * t1015 - t781 * t1030) * t987;
t752 = -t785 * t956 + (t794 * t1016 - t779 * t1031) * t988;
t1 = [t755 * t955 + t756 * t953 + t757 * t951 + m(4) + (-t764 * t964 - t765 * t962 - t766 * t960 + (t772 * t1023 + t771 * t1025 + t770 * t1027) * t929) * t931, -t755 * t956 - t756 * t954 - t757 * t952 + (-t764 * t965 - t765 * t963 - t766 * t961 + (t772 * t1024 + t771 * t1026 + t770 * t1028) * t929) * t931, -t755 * t1004 - t756 * t1002 - t757 * t1000 + (t766 * t944 + t765 * t945 + t764 * t946 + (t772 * t1017 + t771 * t1018 + t770 * t1019) * t929) * t931; t752 * t955 + t753 * t953 + t754 * t951 + (-t761 * t964 - t762 * t962 - t763 * t960 + (t769 * t1023 + t768 * t1025 + t767 * t1027) * t929) * t931, -t752 * t956 - t753 * t954 - t754 * t952 + m(4) + (-t761 * t965 - t762 * t963 - t763 * t961 + (t769 * t1024 + t768 * t1026 + t767 * t1028) * t929) * t931, -t752 * t1004 - t753 * t1002 - t754 * t1000 + (t763 * t944 + t762 * t945 + t761 * t946 + (t769 * t1017 + t768 * t1018 + t767 * t1019) * t929) * t931; t758 * t955 + t759 * t953 + t760 * t951 + (-t773 * t964 - t774 * t962 - t775 * t960 + (t778 * t1023 + t777 * t1025 + t776 * t1027) * t929) * t931, -t758 * t956 - t759 * t954 - t760 * t952 + (-t773 * t965 - t774 * t963 - t775 * t961 + (t778 * t1024 + t777 * t1026 + t776 * t1028) * t929) * t931, -t758 * t1004 - t759 * t1002 - t760 * t1000 + m(4) + (t775 * t944 + t774 * t945 + t773 * t946 + (t778 * t1017 + t777 * t1018 + t776 * t1019) * t929) * t931;];
MX  = t1;
