% Calculate vector of centrifugal and coriolis load on the joints for
% P4PRRRR1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% taucX [4x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:58:30
% EndTime: 2020-08-07 10:58:36
% DurationCPUTime: 6.88s
% Computational Cost: add. (6552->372), mult. (14238->785), div. (5748->13), fcn. (12372->34), ass. (0->331)
t872 = cos(qJ(3,4));
t854 = 0.1e1 / t872;
t914 = 0.1e1 / pkin(2);
t1035 = t854 * t914;
t897 = xP(4);
t850 = sin(t897);
t851 = cos(t897);
t904 = koppelP(4,2);
t908 = koppelP(4,1);
t814 = t850 * t908 + t851 * t904;
t818 = -t850 * t904 + t851 * t908;
t874 = legFrame(4,2);
t842 = sin(t874);
t846 = cos(t874);
t893 = xDP(4);
t895 = xDP(2);
t896 = xDP(1);
t947 = (t818 * t893 + t895) * t842 - (-t814 * t893 + t896) * t846;
t766 = t947 * t1035;
t762 = t766 ^ 2;
t1057 = t762 * t854;
t870 = sin(qJ(3,4));
t1007 = t870 * t1057;
t1078 = rSges(3,3) * m(3);
t840 = rSges(3,1) * t1078 - Icges(3,5);
t1075 = t840 / 0.4e1;
t839 = rSges(3,2) * t1078 - Icges(3,6);
t1076 = -t839 / 0.4e1;
t871 = sin(qJ(2,4));
t1025 = t871 * t872;
t852 = 0.1e1 / t871;
t894 = xDP(3);
t855 = 0.1e1 / t872 ^ 2;
t873 = cos(qJ(2,4));
t989 = t855 * t870 * t873;
t741 = (-t854 * t894 - t947 * t989) * t914 * t852;
t1064 = t741 * t873;
t690 = ((-t766 * t870 * t871 + t872 * t1064) * t854 * t741 + (-t741 * t870 * t1025 + t766 * t873) * t855 * t766) * t852;
t740 = t741 ^ 2;
t710 = (-t740 * t872 - t1057) * t852 * pkin(2);
t912 = rSges(3,2) ^ 2;
t913 = rSges(3,1) ^ 2;
t832 = (-t912 + t913) * m(3) - Icges(3,1) + Icges(3,2);
t1077 = t832 / 0.2e1;
t841 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t900 = 0.2e1 * qJ(3,4);
t1021 = (t912 + t913);
t915 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + ((2 * rSges(3,3) ^ 2 + t1021) * m(3)) / 0.2e1;
t770 = cos(t900) * t1077 - t841 * sin(t900) + t915;
t1092 = (rSges(3,1) * t872 - rSges(3,2) * t870) * m(3);
t838 = m(2) * rSges(2,2) - t1078;
t892 = m(2) * rSges(2,1);
t790 = (t892 + t1092) * t873 - t871 * t838;
t804 = -t839 * t872 - t840 * t870;
t853 = t872 ^ 2;
t995 = t832 * t870 * t872;
t1084 = -t804 * t1007 + t690 * t770 + t710 * t790 + 0.4e1 * ((t872 * t1075 + t870 * t1076) * t766 + (t995 / 0.2e1 + (t853 - 0.1e1 / 0.2e1) * t841) * t741) * t766;
t884 = cos(qJ(3,3));
t861 = 0.1e1 / t884;
t1028 = t861 * t914;
t905 = koppelP(3,2);
t909 = koppelP(3,1);
t815 = t850 * t909 + t851 * t905;
t819 = -t850 * t905 + t851 * t909;
t875 = legFrame(3,2);
t843 = sin(t875);
t847 = cos(t875);
t946 = (t819 * t893 + t895) * t843 - (-t815 * t893 + t896) * t847;
t767 = t946 * t1028;
t763 = t767 ^ 2;
t1056 = t763 * t861;
t878 = sin(qJ(3,3));
t1006 = t878 * t1056;
t879 = sin(qJ(2,3));
t1024 = t879 * t884;
t857 = 0.1e1 / t879;
t862 = 0.1e1 / t884 ^ 2;
t885 = cos(qJ(2,3));
t982 = t862 * t878 * t885;
t751 = (-t861 * t894 - t946 * t982) * t914 * t857;
t1060 = t751 * t885;
t691 = ((-t767 * t878 * t879 + t884 * t1060) * t861 * t751 + (-t751 * t878 * t1024 + t767 * t885) * t862 * t767) * t857;
t742 = t751 ^ 2;
t711 = (-t742 * t884 - t1056) * t857 * pkin(2);
t901 = 0.2e1 * qJ(3,3);
t771 = cos(t901) * t1077 - t841 * sin(t901) + t915;
t1093 = (rSges(3,1) * t884 - rSges(3,2) * t878) * m(3);
t791 = (t892 + t1093) * t885 - t879 * t838;
t805 = -t839 * t884 - t840 * t878;
t860 = t884 ^ 2;
t994 = t832 * t878 * t884;
t1085 = -t805 * t1006 + t691 * t771 + t711 * t791 + 0.4e1 * ((t884 * t1075 + t878 * t1076) * t767 + (t994 / 0.2e1 + (t860 - 0.1e1 / 0.2e1) * t841) * t751) * t767;
t886 = cos(qJ(3,2));
t864 = 0.1e1 / t886;
t1027 = t864 * t914;
t906 = koppelP(2,2);
t910 = koppelP(2,1);
t816 = t850 * t910 + t851 * t906;
t820 = -t850 * t906 + t851 * t910;
t876 = legFrame(2,2);
t844 = sin(t876);
t848 = cos(t876);
t945 = (t820 * t893 + t895) * t844 - (-t816 * t893 + t896) * t848;
t768 = t945 * t1027;
t764 = t768 ^ 2;
t1055 = t764 * t864;
t880 = sin(qJ(3,2));
t1005 = t880 * t1055;
t881 = sin(qJ(2,2));
t1023 = t881 * t886;
t858 = 0.1e1 / t881;
t865 = 0.1e1 / t886 ^ 2;
t887 = cos(qJ(2,2));
t981 = t865 * t880 * t887;
t752 = (-t864 * t894 - t945 * t981) * t914 * t858;
t1059 = t752 * t887;
t692 = ((-t768 * t880 * t881 + t886 * t1059) * t864 * t752 + (-t752 * t880 * t1023 + t768 * t887) * t865 * t768) * t858;
t743 = t752 ^ 2;
t712 = (-t743 * t886 - t1055) * t858 * pkin(2);
t902 = 0.2e1 * qJ(3,2);
t772 = cos(t902) * t1077 - t841 * sin(t902) + t915;
t1094 = (rSges(3,1) * t886 - rSges(3,2) * t880) * m(3);
t792 = (t892 + t1094) * t887 - t881 * t838;
t806 = -t839 * t886 - t840 * t880;
t863 = t886 ^ 2;
t993 = t832 * t880 * t886;
t1086 = -t806 * t1005 + t692 * t772 + t712 * t792 + 0.4e1 * ((t886 * t1075 + t880 * t1076) * t768 + (t993 / 0.2e1 + (t863 - 0.1e1 / 0.2e1) * t841) * t752) * t768;
t888 = cos(qJ(3,1));
t867 = 0.1e1 / t888;
t1026 = t867 * t914;
t907 = koppelP(1,2);
t911 = koppelP(1,1);
t817 = t850 * t911 + t851 * t907;
t821 = -t850 * t907 + t851 * t911;
t877 = legFrame(1,2);
t845 = sin(t877);
t849 = cos(t877);
t944 = (t821 * t893 + t895) * t845 - (-t817 * t893 + t896) * t849;
t769 = t944 * t1026;
t765 = t769 ^ 2;
t1054 = t765 * t867;
t882 = sin(qJ(3,1));
t1004 = t882 * t1054;
t883 = sin(qJ(2,1));
t1022 = t883 * t888;
t859 = 0.1e1 / t883;
t868 = 0.1e1 / t888 ^ 2;
t889 = cos(qJ(2,1));
t980 = t868 * t882 * t889;
t753 = (-t867 * t894 - t944 * t980) * t914 * t859;
t1058 = t753 * t889;
t693 = ((-t769 * t882 * t883 + t888 * t1058) * t867 * t753 + (-t753 * t882 * t1022 + t769 * t889) * t868 * t769) * t859;
t744 = t753 ^ 2;
t713 = (-t744 * t888 - t1054) * t859 * pkin(2);
t903 = 0.2e1 * qJ(3,1);
t773 = cos(t903) * t1077 - t841 * sin(t903) + t915;
t1095 = (rSges(3,1) * t888 - rSges(3,2) * t882) * m(3);
t793 = (t892 + t1095) * t889 - t883 * t838;
t807 = -t839 * t888 - t840 * t882;
t866 = t888 ^ 2;
t992 = t832 * t882 * t888;
t1087 = -t807 * t1004 + t693 * t773 + t713 * t793 + 0.4e1 * ((t888 * t1075 + t882 * t1076) * t769 + (t992 / 0.2e1 + (t866 - 0.1e1 / 0.2e1) * t841) * t753) * t769;
t822 = rSges(3,1) * t870 + rSges(3,2) * t872;
t1074 = m(3) * t822;
t1018 = t871 * t1074;
t1079 = 0.2e1 * t841;
t837 = t1021 * m(3) + Icges(3,3);
t1103 = t837 * t1007 + t710 * t1018 - t690 * t804 + t740 * (t853 * t1079 - t841 + t995);
t823 = rSges(3,1) * t878 + rSges(3,2) * t884;
t1073 = m(3) * t823;
t1016 = t879 * t1073;
t1102 = t837 * t1006 + t711 * t1016 - t691 * t805 + t742 * (t860 * t1079 - t841 + t994);
t824 = rSges(3,1) * t880 + rSges(3,2) * t886;
t1072 = m(3) * t824;
t1014 = t881 * t1072;
t1101 = t837 * t1005 + t712 * t1014 - t692 * t806 + t743 * (t863 * t1079 - t841 + t993);
t825 = rSges(3,1) * t882 + rSges(3,2) * t888;
t1071 = m(3) * t825;
t1012 = t883 * t1071;
t1100 = t837 * t1004 + t713 * t1012 - t693 * t807 + t744 * (t866 * t1079 - t841 + t992);
t1020 = 2 * m(3);
t856 = m(1) + m(2) + m(3);
t979 = t854 * t1018;
t1091 = -t762 * t870 * t979 - t690 * t790 - t710 * t856 + (-t740 * t892 - (t740 + t762) * t1092) * t871 - (t822 * t766 * t1020 + t741 * t838) * t1064;
t978 = t861 * t1016;
t1090 = -t763 * t878 * t978 - t691 * t791 - t711 * t856 + (-t742 * t892 - (t742 + t763) * t1093) * t879 - (t823 * t767 * t1020 + t751 * t838) * t1060;
t977 = t864 * t1014;
t1089 = -t764 * t880 * t977 - t692 * t792 - t712 * t856 + (-t743 * t892 - (t743 + t764) * t1094) * t881 - (t824 * t768 * t1020 + t752 * t838) * t1059;
t976 = t867 * t1012;
t1088 = -t765 * t882 * t976 - t693 * t793 - t713 * t856 + (-t744 * t892 - (t744 + t765) * t1095) * t883 - (t825 * t769 * t1020 + t753 * t838) * t1058;
t964 = t914 * t980;
t948 = t859 * t964;
t1083 = t1100 * t1026 + t1087 * t948;
t965 = t914 * t981;
t949 = t858 * t965;
t1082 = t1101 * t1027 + t1086 * t949;
t966 = t914 * t982;
t950 = t857 * t966;
t1081 = t1102 * t1028 + t1085 * t950;
t970 = t914 * t989;
t951 = t852 * t970;
t1080 = t1103 * t1035 + t1084 * t951;
t869 = t893 ^ 2;
t1070 = m(4) * t869;
t802 = t842 * t1025 - t846 * t870;
t1053 = t802 * t854;
t803 = t846 * t1025 + t842 * t870;
t1052 = t803 * t854;
t808 = t843 * t1024 - t847 * t878;
t1051 = t808 * t861;
t809 = t844 * t1023 - t848 * t880;
t1050 = t809 * t864;
t810 = t845 * t1022 - t849 * t882;
t1049 = t810 * t867;
t811 = t847 * t1024 + t843 * t878;
t1048 = t811 * t861;
t812 = t848 * t1023 + t844 * t880;
t1047 = t812 * t864;
t813 = t849 * t1022 + t845 * t882;
t1046 = t813 * t867;
t1045 = t842 * t914;
t1044 = t843 * t914;
t1043 = t844 * t914;
t1042 = t845 * t914;
t1041 = t846 * t914;
t1040 = t847 * t914;
t1039 = t848 * t914;
t1038 = t849 * t914;
t1037 = t852 * t854;
t1034 = t857 * t861;
t1032 = t858 * t864;
t1030 = t859 * t867;
t1019 = t854 * t1074;
t1017 = t861 * t1073;
t1015 = t864 * t1072;
t1013 = t867 * t1071;
t1003 = t802 * t1037;
t1002 = t803 * t1037;
t1001 = t808 * t1034;
t1000 = t809 * t1032;
t999 = t810 * t1030;
t998 = t811 * t1034;
t997 = t812 * t1032;
t996 = t813 * t1030;
t991 = t856 * t1037;
t990 = t852 * t1035;
t988 = t856 * t1034;
t987 = t856 * t1032;
t986 = t856 * t1030;
t985 = t857 * t1028;
t984 = t858 * t1027;
t983 = t859 * t1026;
t971 = t852 * t989;
t969 = t857 * t982;
t968 = t858 * t981;
t967 = t859 * t980;
t959 = t842 * t970;
t958 = t843 * t966;
t957 = t844 * t965;
t956 = t845 * t964;
t955 = t846 * t970;
t954 = t847 * t966;
t953 = t848 * t965;
t952 = t849 * t964;
t943 = t814 * t846 + t818 * t842;
t942 = t815 * t847 + t819 * t843;
t941 = t816 * t848 + t820 * t844;
t940 = t817 * t849 + t821 * t845;
t931 = (-t814 * t842 + t818 * t846) * t1035;
t930 = (-t815 * t843 + t819 * t847) * t1028;
t929 = (-t816 * t844 + t820 * t848) * t1027;
t928 = (-t817 * t845 + t821 * t849) * t1026;
t927 = t770 * t971 - t804 * t854;
t926 = t771 * t969 - t805 * t861;
t925 = t772 * t968 - t806 * t864;
t924 = t773 * t967 - t807 * t867;
t923 = t804 * t971 - t837 * t854;
t922 = t805 * t969 - t837 * t861;
t921 = t806 * t968 - t837 * t864;
t920 = t807 * t967 - t837 * t867;
t919 = t790 * t971 + t979;
t918 = t791 * t969 + t978;
t917 = t792 * t968 + t977;
t916 = t793 * t967 + t976;
t899 = rSges(4,1);
t898 = rSges(4,2);
t781 = (-t793 * t1026 + t856 * t889) * t859;
t780 = (-t792 * t1027 + t856 * t887) * t858;
t779 = (-t791 * t1028 + t856 * t885) * t857;
t778 = t940 * t1026;
t777 = t941 * t1027;
t776 = t942 * t1028;
t775 = t943 * t1035;
t774 = (-t790 * t1035 + t856 * t873) * t852;
t761 = t940 * t948;
t760 = t941 * t949;
t759 = t942 * t950;
t758 = t943 * t951;
t757 = (-t810 * t817 + t813 * t821) * t1030;
t756 = (-t809 * t816 + t812 * t820) * t1032;
t755 = (-t808 * t815 + t811 * t819) * t1034;
t754 = (-t802 * t814 + t803 * t818) * t1037;
t737 = (-t773 * t1026 + t793 * t889) * t859;
t736 = (-t772 * t1027 + t792 * t887) * t858;
t735 = (-t771 * t1028 + t791 * t885) * t857;
t734 = (-t770 * t1035 + t790 * t873) * t852;
t733 = -t916 * t1042 + t813 * t986;
t732 = -t917 * t1043 + t812 * t987;
t731 = -t918 * t1044 + t811 * t988;
t730 = t916 * t1038 + t810 * t986;
t729 = t917 * t1039 + t809 * t987;
t728 = t918 * t1040 + t808 * t988;
t727 = -t919 * t1045 + t803 * t991;
t726 = t919 * t1041 + t802 * t991;
t722 = -t924 * t1042 + t793 * t996;
t721 = -t925 * t1043 + t792 * t997;
t720 = -t926 * t1044 + t791 * t998;
t719 = t924 * t1038 + t793 * t999;
t718 = t792 * t1000 + t925 * t1039;
t717 = t791 * t1001 + t926 * t1040;
t715 = t790 * t1002 - t927 * t1045;
t714 = t790 * t1003 + t927 * t1041;
t706 = -t778 * t1012 + t757 * t856 - t761 * t793;
t705 = -t777 * t1014 + t756 * t856 - t760 * t792;
t704 = -t776 * t1016 + t755 * t856 - t759 * t791;
t702 = -t775 * t1018 + t754 * t856 - t758 * t790;
t701 = t757 * t793 - t761 * t773 + t778 * t807;
t700 = t756 * t792 - t760 * t772 + t777 * t806;
t699 = t755 * t791 - t759 * t771 + t776 * t805;
t698 = t754 * t790 - t758 * t770 + t775 * t804;
t1 = [(t850 * t898 - t851 * t899) * t1070 + t1088 * t999 + t1091 * t1003 + t1090 * t1001 + t1089 * t1000 + ((-t802 * t1019 + t923 * t1041) * t931 + (-(t726 * t1053 + t714 * t955) * t818 - (t726 * t1052 - t714 * t959) * t814) * t852 + (-t808 * t1017 + t922 * t1040) * t930 + (-(t728 * t1051 + t717 * t954) * t819 - (t728 * t1048 - t717 * t958) * t815) * t857 + (-t810 * t1013 + t920 * t1038) * t928 + (-(t730 * t1049 + t719 * t952) * t821 - (t730 * t1046 - t719 * t956) * t817) * t859 + (-t809 * t1015 + t921 * t1039) * t929 + (-(t729 * t1050 + t718 * t953) * t820 - (t729 * t1047 - t718 * t957) * t816) * t858) * t869 - t1083 * t849 - t1082 * t848 - t1081 * t847 - t1080 * t846; -(t850 * t899 + t851 * t898) * t1070 + t1090 * t998 + t1089 * t997 + t1088 * t996 + t1091 * t1002 + ((-t803 * t1019 - t923 * t1045) * t931 + (-(t727 * t1053 + t715 * t955) * t818 - (t727 * t1052 - t715 * t959) * t814) * t852 + (-t811 * t1017 - t922 * t1044) * t930 + (-(t731 * t1051 + t720 * t954) * t819 - (t731 * t1048 - t720 * t958) * t815) * t857 + (-t813 * t1013 - t920 * t1042) * t928 + (-(t733 * t1049 + t722 * t952) * t821 - (t733 * t1046 - t722 * t956) * t817) * t859 + (-t812 * t1015 - t921 * t1043) * t929 + (-(t732 * t1050 + t721 * t953) * t820 - (t732 * t1047 - t721 * t957) * t816) * t858) * t869 + t1083 * t845 + t1082 * t844 + t1081 * t843 + t1080 * t842; t1084 * t990 + t1085 * t985 + t1086 * t984 + t1087 * t983 + t1091 * t852 * t873 + t1090 * t857 * t885 + t1089 * t858 * t887 + t1088 * t859 * t889 + ((-t889 * t1071 - t807 * t983) * t928 + (-(t781 * t1049 + t737 * t952) * t821 - (t781 * t1046 - t737 * t956) * t817) * t859 + (-t887 * t1072 - t806 * t984) * t929 + (-(t780 * t1050 + t736 * t953) * t820 - (t780 * t1047 - t736 * t957) * t816) * t858 + (-t885 * t1073 - t805 * t985) * t930 + (-(t779 * t1051 + t735 * t954) * t819 - (t779 * t1048 - t735 * t958) * t815) * t857 + (-t873 * t1074 - t804 * t990) * t931 + (-(t774 * t1053 + t734 * t955) * t818 - (t774 * t1052 - t734 * t959) * t814) * t852) * t869; t1100 * t778 + t1101 * t777 + t1102 * t776 + t1103 * t775 + t1087 * t761 + t1086 * t760 + t1085 * t759 + t1084 * t758 + t1088 * t757 + t1089 * t756 + t1090 * t755 + t1091 * t754 + ((-t755 * t1016 - t759 * t805 + t776 * t837) * t930 + (-(t704 * t1051 + t699 * t954) * t819 - (t704 * t1048 - t699 * t958) * t815) * t857 + (-t757 * t1012 - t761 * t807 + t778 * t837) * t928 + (-(t706 * t1049 + t701 * t952) * t821 - (t706 * t1046 - t701 * t956) * t817) * t859 + (-t754 * t1018 - t758 * t804 + t775 * t837) * t931 + (-(t702 * t1053 + t698 * t955) * t818 - (t702 * t1052 - t698 * t959) * t814) * t852 + (-t756 * t1014 - t760 * t806 + t777 * t837) * t929 + (-(t705 * t1050 + t700 * t953) * t820 - (t705 * t1047 - t700 * t957) * t816) * t858) * t869;];
taucX  = t1;
