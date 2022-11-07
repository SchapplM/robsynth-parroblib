% Calculate vector of centrifugal and Coriolis load for parallel robot
% P3RRPRR8V2G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2022-11-07 13:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:05:44
% EndTime: 2022-11-07 13:05:49
% DurationCPUTime: 5.20s
% Computational Cost: add. (26829->449), mult. (46416->692), div. (4257->12), fcn. (33177->50), ass. (0->332)
t907 = rSges(2,3) + pkin(5);
t1087 = 2 * pkin(1);
t987 = 2 * pkin(3);
t1086 = pkin(1) / 0.2e1;
t1085 = pkin(2) / 0.2e1;
t1084 = 2 * m(3);
t875 = qJ(2,3) + pkin(7);
t844 = cos(t875);
t821 = pkin(3) * t844;
t901 = cos(qJ(2,3));
t857 = t901 * pkin(2);
t994 = t857 + t821;
t793 = 0.1e1 / t994;
t892 = legFrame(3,2);
t850 = sin(t892);
t853 = cos(t892);
t914 = xDP(2);
t915 = xDP(1);
t765 = (t850 * t915 + t853 * t914) * t793;
t762 = t765 ^ 2;
t876 = qJ(2,2) + pkin(7);
t845 = cos(t876);
t822 = pkin(3) * t845;
t903 = cos(qJ(2,2));
t858 = t903 * pkin(2);
t993 = t858 + t822;
t794 = 0.1e1 / t993;
t893 = legFrame(2,2);
t851 = sin(t893);
t854 = cos(t893);
t766 = (t851 * t915 + t854 * t914) * t794;
t763 = t766 ^ 2;
t877 = qJ(2,1) + pkin(7);
t846 = cos(t877);
t823 = pkin(3) * t846;
t905 = cos(qJ(2,1));
t859 = t905 * pkin(2);
t992 = t859 + t823;
t795 = 0.1e1 / t992;
t894 = legFrame(1,2);
t852 = sin(t894);
t855 = cos(t894);
t767 = (t852 * t915 + t855 * t914) * t795;
t764 = t767 ^ 2;
t890 = cos(pkin(7));
t1041 = t890 * pkin(3);
t824 = pkin(2) + t1041;
t889 = sin(pkin(7));
t895 = sin(qJ(2,3));
t1006 = t889 * t895;
t968 = pkin(3) * t1006;
t942 = -t824 * t901 + t968;
t787 = 0.1e1 / t942;
t1038 = pkin(5) + qJ(3,3);
t867 = pkin(6) + t1038;
t847 = 0.1e1 / t867;
t1033 = t787 * t847;
t896 = sin(qJ(1,3));
t1014 = t850 * t896;
t1015 = t850 * t824;
t1024 = t824 * t853;
t902 = cos(qJ(1,3));
t958 = pkin(1) * t896 - t902 * t867;
t871 = t890 ^ 2;
t965 = pkin(3) * (t871 - 0.1e1);
t1045 = (t958 * t1006 + t896 * t965) * pkin(3);
t1077 = 0.2e1 * t901 ^ 2;
t1042 = t889 * pkin(3);
t833 = pkin(1) * t1042;
t832 = pkin(2) * t1041;
t1078 = 0.2e1 * t832;
t931 = pkin(3) ^ 2;
t932 = pkin(2) ^ 2;
t937 = 0.2e1 * t871 * t931 + t1078 - t931 + t932;
t775 = t937 * t895 + t833;
t781 = -0.2e1 * t896 * t968 + t958;
t891 = t932 / 0.2e1;
t998 = t832 + t891;
t790 = (t871 - 0.1e1 / 0.2e1) * t931 + t998;
t800 = pkin(1) * t895 - t1042;
t974 = t853 * t1042;
t741 = (-t790 * t1014 + t824 * t974) * t1077 + (-t781 * t1015 + t853 * t775) * t901 + t850 * t1045 + t800 * t1024;
t735 = t741 * t914 * t1033;
t1009 = t853 * t896;
t971 = t850 * t1042;
t742 = (t790 * t1009 + t824 * t971) * t1077 + (t781 * t1024 + t775 * t850) * t901 - t853 * t1045 + t800 * t1015;
t736 = t742 * t915 * t1033;
t1036 = t857 + pkin(1);
t778 = t896 * t867 + (t1036 + t821) * t902;
t913 = xDP(3);
t768 = t778 * t847 * t913;
t729 = -t736 - t735 + t768;
t1082 = 0.2e1 * t729;
t897 = sin(qJ(2,2));
t1005 = t889 * t897;
t967 = pkin(3) * t1005;
t941 = -t824 * t903 + t967;
t788 = 0.1e1 / t941;
t1039 = pkin(5) + qJ(3,2);
t868 = pkin(6) + t1039;
t848 = 0.1e1 / t868;
t1032 = t788 * t848;
t898 = sin(qJ(1,2));
t1012 = t851 * t898;
t1013 = t851 * t824;
t1023 = t824 * t854;
t904 = cos(qJ(1,2));
t957 = pkin(1) * t898 - t904 * t868;
t1044 = (t957 * t1005 + t898 * t965) * pkin(3);
t1076 = 0.2e1 * t903 ^ 2;
t776 = t937 * t897 + t833;
t782 = -0.2e1 * t898 * t967 + t957;
t801 = pkin(1) * t897 - t1042;
t973 = t854 * t1042;
t743 = (-t790 * t1012 + t824 * t973) * t1076 + (-t782 * t1013 + t854 * t776) * t903 + t851 * t1044 + t801 * t1023;
t737 = t743 * t914 * t1032;
t1008 = t854 * t898;
t970 = t851 * t1042;
t744 = (t790 * t1008 + t824 * t970) * t1076 + (t782 * t1023 + t776 * t851) * t903 - t854 * t1044 + t801 * t1013;
t738 = t744 * t915 * t1032;
t1035 = t858 + pkin(1);
t779 = t898 * t868 + (t1035 + t822) * t904;
t769 = t779 * t848 * t913;
t730 = -t738 - t737 + t769;
t1081 = 0.2e1 * t730;
t899 = sin(qJ(2,1));
t1004 = t889 * t899;
t966 = pkin(3) * t1004;
t940 = -t824 * t905 + t966;
t789 = 0.1e1 / t940;
t1037 = qJ(3,1) + pkin(5);
t869 = pkin(6) + t1037;
t849 = 0.1e1 / t869;
t1031 = t789 * t849;
t900 = sin(qJ(1,1));
t1010 = t852 * t900;
t1011 = t852 * t824;
t1022 = t824 * t855;
t906 = cos(qJ(1,1));
t956 = pkin(1) * t900 - t906 * t869;
t1043 = (t956 * t1004 + t900 * t965) * pkin(3);
t1075 = 0.2e1 * t905 ^ 2;
t777 = t937 * t899 + t833;
t783 = -0.2e1 * t900 * t966 + t956;
t802 = pkin(1) * t899 - t1042;
t972 = t855 * t1042;
t745 = (-t790 * t1010 + t824 * t972) * t1075 + (-t783 * t1011 + t855 * t777) * t905 + t852 * t1043 + t802 * t1022;
t739 = t745 * t914 * t1031;
t1007 = t855 * t900;
t969 = t852 * t1042;
t746 = (t790 * t1007 + t824 * t969) * t1075 + (t783 * t1022 + t777 * t852) * t905 - t855 * t1043 + t802 * t1011;
t740 = t746 * t915 * t1031;
t1034 = t859 + pkin(1);
t780 = t900 * t869 + (t1034 + t823) * t906;
t770 = t780 * t849 * t913;
t731 = -t740 - t739 + t770;
t1080 = 0.2e1 * t731;
t1079 = -0.4e1 * pkin(1) * (t931 / 0.2e1 + t998);
t1074 = -0.4e1 * t901;
t1073 = -0.4e1 * t903;
t1072 = -0.4e1 * t905;
t1071 = -4 * pkin(5) - 4 * pkin(6);
t1070 = pkin(2) * m(3);
t1069 = m(2) * rSges(2,1);
t1068 = m(2) * rSges(2,2);
t1067 = m(3) * rSges(3,2);
t1066 = t793 / 0.2e1;
t1065 = t794 / 0.2e1;
t1064 = t795 / 0.2e1;
t926 = rSges(2,2) ^ 2;
t928 = rSges(2,1) ^ 2;
t796 = t932 * m(3) + (-t926 + t928) * m(2) + Icges(2,2) - Icges(2,1);
t1063 = -t796 / 0.2e1;
t925 = rSges(3,2) ^ 2;
t927 = (rSges(3,1) ^ 2);
t805 = (-t925 + t927) * m(3) - Icges(3,1) + Icges(3,2);
t1062 = -t805 / 0.2e1;
t829 = -rSges(3,1) * t1067 + Icges(3,4);
t1061 = -t829 / 0.2e1;
t830 = rSges(2,1) * t1068 - Icges(2,4);
t1060 = t830 / 0.2e1;
t1059 = m(2) * t907;
t864 = (rSges(3,3) + t1038);
t1058 = m(3) * t864;
t865 = (rSges(3,3) + t1039);
t1057 = m(3) * t865;
t866 = (rSges(3,3) + t1037);
t1056 = m(3) * t866;
t1055 = pkin(1) * t844;
t1054 = pkin(1) * t845;
t1053 = pkin(1) * t846;
t922 = 0.2e1 * qJ(2,3);
t872 = pkin(7) + t922;
t835 = sin(t872);
t1052 = pkin(2) * t835;
t923 = 0.2e1 * qJ(2,2);
t873 = pkin(7) + t923;
t836 = sin(t873);
t1051 = pkin(2) * t836;
t924 = 0.2e1 * qJ(2,1);
t874 = pkin(7) + t924;
t837 = sin(t874);
t1050 = pkin(2) * t837;
t1049 = pkin(2) * t895;
t1048 = pkin(2) * t897;
t1047 = pkin(2) * t899;
t1046 = pkin(3) * t932;
t1040 = t931 * pkin(2);
t879 = sin(t922);
t1030 = t796 * t879;
t880 = sin(t923);
t1029 = t796 * t880;
t881 = sin(t924);
t1028 = t796 * t881;
t825 = 0.2e1 * t875;
t813 = sin(t825);
t1027 = t805 * t813;
t826 = 0.2e1 * t876;
t814 = sin(t826);
t1026 = t805 * t814;
t827 = 0.2e1 * t877;
t815 = sin(t827);
t1025 = t805 * t815;
t816 = cos(t825);
t1021 = t829 * t816;
t817 = cos(t826);
t1020 = t829 * t817;
t818 = cos(t827);
t1019 = t829 * t818;
t882 = cos(t922);
t1018 = t830 * t882;
t883 = cos(t923);
t1017 = t830 * t883;
t884 = cos(t924);
t1016 = t830 * t884;
t1003 = pkin(2) * t987;
t1002 = rSges(2,1) * t1059 - Icges(2,5);
t756 = (-t824 * t1014 + t974) * t901 + (t896 * t971 + t1024) * t895;
t759 = (t824 * t1009 + t971) * t901 + t895 * (-t896 * t974 + t1015);
t750 = (t902 * t913 - (t756 * t914 + t759 * t915) * t787) * t847;
t747 = pkin(1) * t750;
t723 = t747 + t736 / 0.2e1 + t735 / 0.2e1 - t768 / 0.2e1;
t838 = sin(t875);
t797 = pkin(3) * t838 + t1049;
t831 = pkin(2) * t891 + t1040;
t921 = 0.2e1 * pkin(7);
t929 = pkin(5) ^ 2;
t933 = pkin(1) ^ 2;
t990 = t931 + t932;
t939 = -(2 * t929) - (2 * t933) - t990 + ((-4 * pkin(5) - 2 * pkin(6)) * pkin(6));
t959 = -(2 * pkin(3) * t931) - 0.4e1 * t1046;
t981 = -0.2e1 * t1040;
t982 = -0.2e1 * t1046;
t841 = cos(t872);
t997 = -t841 - t890;
t708 = ((cos(-pkin(7) + qJ(2,3)) * t982 + cos(t921 + qJ(2,3)) * t981 + t959 * t844 + t831 * t1074 + t1079) * t762 * t1066 + ((pkin(1) + t994) * t729 - 0.2e1 * t723 * t821 + t1082 * t1086 + (-t931 * t816 - t932 * t882 + t939 + ((t1071 - 2 * qJ(3,3)) * qJ(3,3))) * t750 / 0.2e1 + (t997 * t750 * t987 + t723 * t1074) * t1085 + (t797 + (t835 * t1003 + t813 * t931 + t879 * t932) * t1066) * t765 * t867) * t750) * t847;
t803 = t1078 + t990;
t720 = (-t762 * t803 / (t857 + (t890 * t901 - t1006) * pkin(3)) + (t942 * t750 + t1082 - t747) * t750) * t847;
t807 = rSges(3,2) * t1058 - Icges(3,6);
t810 = rSges(3,1) * t1058 - Icges(3,5);
t946 = -t907 * t1069 + Icges(2,5);
t960 = t907 * t1068 - Icges(2,6);
t978 = pkin(2) * t1058;
t753 = -(t807 * t889 - t810 * t890 + t946 - t978) * t895 + (t807 * t890 + t810 * t889 + t960) * t901;
t784 = -t844 * rSges(3,1) + t838 * rSges(3,2) - t1036;
t819 = rSges(2,2) * t1059 - Icges(2,6);
t820 = t889 * pkin(2) * t1067;
t919 = 2 * t933;
t991 = t926 + t928;
t936 = -(t919 + (2 * t929) + ((4 * pkin(5) + 2 * rSges(2,3)) * rSges(2,3)) + t991) * m(2) / 0.2e1 - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t820 - Icges(3,2) / 0.2e1 - Icges(2,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,1) / 0.2e1 - Icges(1,3);
t938 = -t919 / 0.2e1 - t925 / 0.2e1 - t927 / 0.2e1 - t932 / 0.2e1;
t828 = t1069 + t1070;
t945 = t901 * t1068 + t828 * t895;
t963 = t762 * t793 * t797;
t964 = t1068 * t1087;
t979 = 2 * rSges(3,2);
t980 = 2 * rSges(3,1);
t983 = -0.2e1 * pkin(1) * t828;
t986 = t750 * t1084;
t989 = rSges(3,2) * t1087;
t1001 = (t816 * t1062 + t882 * t1063 + t895 * t964 + t901 * t983 + t936) * t720 - t753 * t963 + 0.2e1 * (t879 * t1060 + t813 * t1061) * t720 + (-t784 * t708 + (rSges(3,2) * t1052 + t838 * t989 - (t864 ^ 2) + t938 + (t997 * pkin(2) - 0.2e1 * t1055) * rSges(3,1)) * t720) * m(3) + t729 * t864 * t986 + ((-t810 * t844 + t807 * t838 - (t978 + t1002) * t901 + t819 * t895) * t765 + (-t1030 - t1027 + 0.2e1 * t1021 - 0.2e1 * t1018 - t945 * t1087 + ((-pkin(2) * t841 - t1055) * t979 + (-pkin(1) * t838 - t1052) * t980) * m(3)) * t750) * t765;
t757 = (-t824 * t1012 + t973) * t903 + (t898 * t970 + t1023) * t897;
t760 = (t824 * t1008 + t970) * t903 + t897 * (-t898 * t973 + t1013);
t751 = (t904 * t913 - (t757 * t914 + t760 * t915) * t788) * t848;
t748 = pkin(1) * t751;
t724 = t748 + t738 / 0.2e1 + t737 / 0.2e1 - t769 / 0.2e1;
t839 = sin(t876);
t799 = pkin(3) * t839 + t1048;
t842 = cos(t873);
t996 = -t842 - t890;
t709 = ((cos(qJ(2,2) - pkin(7)) * t982 + cos(t921 + qJ(2,2)) * t981 + t959 * t845 + t831 * t1073 + t1079) * t763 * t1065 + ((pkin(1) + t993) * t730 - 0.2e1 * t724 * t822 + t1081 * t1086 + (-t931 * t817 - t932 * t883 + t939 + ((t1071 - 2 * qJ(3,2)) * qJ(3,2))) * t751 / 0.2e1 + (t996 * t751 * t987 + t724 * t1073) * t1085 + (t799 + (t836 * t1003 + t814 * t931 + t880 * t932) * t1065) * t766 * t868) * t751) * t848;
t721 = (-t763 * t803 / (t858 + (t890 * t903 - t1005) * pkin(3)) + (t941 * t751 + t1081 - t748) * t751) * t848;
t808 = rSges(3,2) * t1057 - Icges(3,6);
t811 = rSges(3,1) * t1057 - Icges(3,5);
t977 = pkin(2) * t1057;
t754 = -(t808 * t889 - t811 * t890 + t946 - t977) * t897 + (t808 * t890 + t811 * t889 + t960) * t903;
t785 = -t845 * rSges(3,1) + t839 * rSges(3,2) - t1035;
t944 = t903 * t1068 + t828 * t897;
t962 = t763 * t794 * t799;
t985 = t751 * t1084;
t1000 = (t817 * t1062 + t883 * t1063 + t897 * t964 + t903 * t983 + t936) * t721 - t754 * t962 + 0.2e1 * (t880 * t1060 + t814 * t1061) * t721 + (-t785 * t709 + (rSges(3,2) * t1051 + t839 * t989 - (t865 ^ 2) + t938 + (t996 * pkin(2) - 0.2e1 * t1054) * rSges(3,1)) * t721) * m(3) + t730 * t865 * t985 + ((-t811 * t845 + t808 * t839 - (t977 + t1002) * t903 + t819 * t897) * t766 + (-t1029 - t1026 + 0.2e1 * t1020 - 0.2e1 * t1017 - t944 * t1087 + ((-pkin(2) * t842 - t1054) * t979 + (-pkin(1) * t839 - t1051) * t980) * m(3)) * t751) * t766;
t758 = (-t824 * t1010 + t972) * t905 + (t900 * t969 + t1022) * t899;
t761 = (t824 * t1007 + t969) * t905 + t899 * (-t900 * t972 + t1011);
t752 = (t906 * t913 - (t758 * t914 + t761 * t915) * t789) * t849;
t749 = pkin(1) * t752;
t725 = t749 + t740 / 0.2e1 + t739 / 0.2e1 - t770 / 0.2e1;
t840 = sin(t877);
t798 = pkin(3) * t840 + t1047;
t843 = cos(t874);
t995 = -t843 - t890;
t710 = ((cos(qJ(2,1) - pkin(7)) * t982 + cos(qJ(2,1) + t921) * t981 + t959 * t846 + t831 * t1072 + t1079) * t764 * t1064 + ((pkin(1) + t992) * t731 - 0.2e1 * t725 * t823 + t1080 * t1086 + (-t931 * t818 - t932 * t884 + t939 + ((t1071 - 2 * qJ(3,1)) * qJ(3,1))) * t752 / 0.2e1 + (t995 * t752 * t987 + t725 * t1072) * t1085 + (t798 + (t837 * t1003 + t815 * t931 + t881 * t932) * t1064) * t767 * t869) * t752) * t849;
t722 = (-t764 * t803 / (t859 + (t890 * t905 - t1004) * pkin(3)) + (t940 * t752 + t1080 - t749) * t752) * t849;
t809 = rSges(3,2) * t1056 - Icges(3,6);
t812 = rSges(3,1) * t1056 - Icges(3,5);
t976 = pkin(2) * t1056;
t755 = -(t809 * t889 - t812 * t890 + t946 - t976) * t899 + (t809 * t890 + t812 * t889 + t960) * t905;
t786 = -t846 * rSges(3,1) + t840 * rSges(3,2) - t1034;
t943 = t905 * t1068 + t828 * t899;
t961 = t764 * t795 * t798;
t984 = t752 * t1084;
t999 = (t818 * t1062 + t884 * t1063 + t899 * t964 + t905 * t983 + t936) * t722 - t755 * t961 + 0.2e1 * (t881 * t1060 + t815 * t1061) * t722 + (-t786 * t710 + (rSges(3,2) * t1050 + t840 * t989 - (t866 ^ 2) + t938 + (t995 * pkin(2) - 0.2e1 * t1053) * rSges(3,1)) * t722) * m(3) + t731 * t866 * t984 + ((-t812 * t846 + t809 * t840 - (t976 + t1002) * t905 + t819 * t899) * t767 + (-t1028 - t1025 + 0.2e1 * t1019 - 0.2e1 * t1016 - t943 * t1087 + ((-pkin(2) * t843 - t1053) * t979 + (-pkin(1) * t840 - t1050) * t980) * m(3)) * t752) * t767;
t952 = rSges(3,1) * t838 + rSges(3,2) * t844;
t955 = (-t750 * t864 / 0.2e1 + (t952 + t1049) * t765) * t986 + (-t720 * t784 - t708) * m(3);
t951 = rSges(3,1) * t839 + rSges(3,2) * t845;
t954 = (-t751 * t865 / 0.2e1 + (t951 + t1048) * t766) * t985 + (-t721 * t785 - t709) * m(3);
t950 = rSges(3,1) * t840 + rSges(3,2) * t846;
t953 = (-t752 * t866 / 0.2e1 + (t950 + t1047) * t767) * t984 + (-t722 * t786 - t710) * m(3);
t774 = 0.2e1 * t820 - t991 * m(2) - Icges(2,3) - Icges(3,3) + (-0.2e1 * t890 * rSges(3,1) * pkin(2) - t925 - t927 - t932) * m(3);
t949 = t793 * (((-0.2e1 * t729 * t1049 + t952 * (t747 + 0.2e1 * t736 + 0.2e1 * t735 - 0.2e1 * t768)) * m(3) + (t1027 / 0.2e1 - t1021 + t1030 / 0.2e1 + t1018 + t945 * pkin(1) + (rSges(3,1) * t835 + rSges(3,2) * t841) * t1070) * t750) * t750 + t753 * t720 - t774 * t963);
t948 = t794 * (((-0.2e1 * t730 * t1048 + t951 * (t748 + 0.2e1 * t738 + 0.2e1 * t737 - 0.2e1 * t769)) * m(3) + (t1026 / 0.2e1 - t1020 + t1029 / 0.2e1 + t1017 + t944 * pkin(1) + (rSges(3,1) * t836 + rSges(3,2) * t842) * t1070) * t751) * t751 + t754 * t721 - t774 * t962);
t947 = t795 * (((-0.2e1 * t731 * t1047 + t950 * (t749 + 0.2e1 * t740 + 0.2e1 * t739 - 0.2e1 * t770)) * m(3) + (t1025 / 0.2e1 - t1019 + t1028 / 0.2e1 + t1016 + t943 * pkin(1) + (rSges(3,1) * t837 + rSges(3,2) * t843) * t1070) * t752) * t752 + t755 * t722 - t774 * t961);
t1 = [t852 * t947 + t851 * t948 + t850 * t949 - (t953 * t746 + t999 * t761) * t1031 - (t1000 * t760 + t954 * t744) * t1032 - (t1001 * t759 + t955 * t742) * t1033; t855 * t947 + t854 * t948 + t853 * t949 - (t953 * t745 + t999 * t758) * t1031 - (t1000 * t757 + t954 * t743) * t1032 - (t1001 * t756 + t955 * t741) * t1033; (t953 * t780 + t999 * t906) * t849 + (t1000 * t904 + t954 * t779) * t848 + (t1001 * t902 + t955 * t778) * t847;];
taucX  = t1;
