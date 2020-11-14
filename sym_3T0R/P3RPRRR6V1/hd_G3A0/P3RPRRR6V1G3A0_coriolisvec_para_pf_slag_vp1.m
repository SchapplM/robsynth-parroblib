% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR6V1G3A0
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:41:21
% EndTime: 2020-08-06 18:41:27
% DurationCPUTime: 6.92s
% Computational Cost: add. (40131->451), mult. (45354->702), div. (6933->16), fcn. (30555->88), ass. (0->313)
t936 = sin(qJ(1,3));
t955 = -pkin(6) - pkin(5);
t1079 = t936 * t955;
t927 = sin(pkin(7));
t1080 = t936 * t927;
t1070 = t955 * t927;
t845 = -pkin(1) + t1070;
t942 = cos(qJ(1,3));
t825 = t845 * t942;
t928 = cos(pkin(7));
t1148 = t825 + pkin(2) * t1080 - (pkin(2) * t942 - t1079) * t928;
t938 = sin(qJ(1,2));
t1075 = t938 * t955;
t1076 = t938 * t927;
t944 = cos(qJ(1,2));
t826 = t845 * t944;
t1147 = t826 + pkin(2) * t1076 - (pkin(2) * t944 - t1075) * t928;
t940 = sin(qJ(1,1));
t1071 = t940 * t955;
t1072 = t940 * t927;
t946 = cos(qJ(1,1));
t827 = t845 * t946;
t1146 = t827 + pkin(2) * t1072 - (pkin(2) * t946 - t1071) * t928;
t937 = sin(qJ(3,2));
t1113 = t937 * pkin(2);
t970 = 0.2e1 * qJ(3,2);
t907 = sin(t970);
t1114 = pkin(3) * t907;
t1129 = 0.2e1 * t955;
t1132 = 0.2e1 * pkin(1);
t893 = pkin(7) + qJ(3,2);
t869 = sin(t893);
t895 = -pkin(7) + qJ(3,2);
t870 = sin(t895);
t900 = qJ(1,2) + pkin(7);
t876 = cos(t900);
t976 = 0.2e1 * pkin(2);
t1098 = (t876 * t1129 + sin(t900) * t976 + t938 * t1132 + (sin(qJ(1,2) - t895) + sin(qJ(1,2) + t893)) * pkin(3)) / (0.2e1 * t1113 + t1114 + (t869 + t870) * pkin(1));
t952 = xDP(3);
t1026 = t952 * t1098;
t943 = cos(qJ(3,2));
t887 = t943 * pkin(3);
t863 = t887 + pkin(2);
t1090 = t863 * t944;
t788 = (-t1075 + t1090) * t928 - t826 - t863 * t1076;
t879 = t928 * pkin(1);
t829 = 0.1e1 / (t879 + t863);
t921 = 0.1e1 / t937;
t1028 = t788 * t829 * t921;
t933 = legFrame(2,2);
t884 = cos(t933);
t954 = xDP(1);
t1086 = t884 * t954;
t881 = sin(t933);
t953 = xDP(2);
t1087 = t881 * t953;
t978 = pkin(3) ^ 2;
t979 = 0.1e1 / pkin(3);
t1145 = (t1026 / 0.4e1 - (t1086 / 0.4e1 - t1087 / 0.4e1) * t1028) * t978 * t979;
t1144 = t953 / 0.2e1;
t1143 = t954 / 0.2e1;
t1101 = t788 * t979;
t1027 = t921 * t1101;
t1142 = (t1086 / 0.6e1 - t1087 / 0.6e1) * t1027;
t1141 = (t1086 / 0.3e1 - t1087 / 0.3e1) * t1027;
t1140 = (t1086 / 0.2e1 - t1087 / 0.2e1) * t1027;
t1118 = pkin(1) * t927;
t1139 = -t955 + t1118;
t856 = t879 + pkin(2);
t1031 = m(3) * t856 / 0.2e1;
t1011 = rSges(3,1) * t1031;
t1127 = m(3) * rSges(3,2);
t1047 = rSges(3,1) * t1127;
t1063 = -t1047 / 0.2e1 + Icges(3,4) / 0.2e1;
t939 = sin(qJ(3,1));
t945 = cos(qJ(3,1));
t840 = t945 * rSges(3,1) - t939 * rSges(3,2);
t1119 = m(3) * t840;
t964 = 0.2e1 * pkin(7);
t890 = cos(t964);
t981 = pkin(1) ^ 2;
t867 = t981 * t890;
t1135 = -t867 - t978 / 0.2e1;
t980 = pkin(2) ^ 2;
t1024 = -0.2e1 * t980 + t1135;
t1007 = -t981 + t1024;
t963 = -0.2e1 * t981;
t1009 = -0.2e1 * t867 - 0.4e1 * t980 + t963 - t978;
t1112 = t955 * pkin(1);
t1038 = t928 * t1112;
t889 = sin(t964);
t1065 = t981 * t889;
t1012 = -0.2e1 * t1038 + t1065;
t1057 = 0.2e1 * pkin(3);
t1023 = t1139 * t1057;
t888 = t945 * pkin(3);
t864 = t888 + t976;
t1040 = t864 * t879;
t1128 = -0.6e1 * t978;
t1043 = t856 * t1128;
t1044 = -0.2e1 * t1139 * t978;
t1045 = pkin(2) * t879;
t1039 = pkin(1) * t1070;
t847 = -0.2e1 * t1039;
t926 = t955 ^ 2;
t958 = 0.3e1 * t980;
t962 = 0.2e1 * t981;
t1053 = -0.4e1 * pkin(3) * (0.6e1 * t1045 + t847 + t962 + t958 + t926 - t1135);
t1056 = 0.4e1 * pkin(3);
t971 = 0.4e1 * qJ(3,1);
t1067 = t978 * cos(t971);
t959 = 0.2e1 * t980;
t816 = 0.4e1 * t1045 + t867 + t959 + t981;
t973 = 0.2e1 * qJ(3,1);
t919 = cos(t973);
t1093 = t816 * t919;
t1116 = pkin(2) * t1139;
t806 = t1012 + 0.2e1 * t1116;
t1095 = t806 * t919;
t1122 = t979 / 0.2e1;
t1130 = -0.2e1 * t978 - 0.2e1 * t816;
t1131 = 0.4e1 * pkin(2);
t1049 = pkin(7) + qJ(3,1);
t1051 = -pkin(7) + qJ(3,1);
t901 = qJ(1,1) + pkin(7);
t877 = cos(t901);
t910 = sin(t973);
t1097 = (t877 * t1129 + sin(t901) * t976 + t940 * t1132 + (sin(qJ(1,1) - t1051) + sin(qJ(1,1) + t1049)) * pkin(3)) / (t939 * t976 + pkin(3) * t910 + (sin(t1049) + sin(t1051)) * pkin(1));
t865 = t888 + pkin(2);
t1089 = t865 * t946;
t790 = (-t1071 + t1089) * t928 - t827 - t865 * t1072;
t830 = 0.1e1 / (t879 + t865);
t934 = legFrame(1,2);
t882 = sin(t934);
t885 = cos(t934);
t922 = 0.1e1 / t939;
t767 = (t952 * t1097 + (t882 * t953 - t885 * t954) * t922 * t830 * t790) * t979;
t1103 = t767 * t939;
t1032 = pkin(3) * t1103;
t854 = t934 + t901;
t855 = -t934 + t901;
t812 = -sin(t854) - sin(t855);
t815 = cos(t855) - cos(t854);
t781 = (t812 * t1143 + t815 * t1144 - t877 * t952) * t830;
t759 = t1139 * t781 - t1032;
t908 = sin(t971);
t972 = 0.3e1 * qJ(3,1);
t909 = sin(t972);
t918 = cos(t972);
t977 = pkin(3) * t978;
t1062 = t926 / 0.2e1 + t981;
t1008 = 0.3e1 / 0.8e1 * t978 + t980 / 0.2e1 + t1062;
t1060 = t926 + t981;
t984 = -0.8e1 * (0.3e1 / 0.4e1 * t978 + t958 + t1060) * t879 - 0.8e1 * (t1008 - t1039) * t976 + 0.8e1 * (-pkin(2) * t890 + t889 * t955) * t981;
t990 = pkin(3) * (-t945 * pkin(2) + t856 * t918);
t1042 = 0.2e1 * t1065;
t995 = -0.4e1 * t1038 + t1042;
t747 = ((t918 * t1044 + (t1139 * t864 + t1012 - t1095) * t1057) * t767 + (t909 * t1043 + t910 * t1053 - t977 * t908 + t984 * t939) * t781) / (t1093 + t1067 / 0.2e1 - 0.2e1 * t1040 + 0.2e1 * t990 + t1007) * t781 * t1122 + (t759 * t1131 + (t910 * t1130 - t978 * t908 + (-t856 * t909 - t939 * t879) * t1056) * t767 + (-0.2e1 * t1095 + (-t918 + t945) * t1023 + t995) * t781) / (t1009 - 0.4e1 * t1040 + t1067 + 0.2e1 * t1093 + 0.4e1 * t990) * t767;
t1083 = t922 * t945;
t764 = t767 ^ 2;
t808 = t847 + t980 + 0.2e1 * t1045 + t1060;
t925 = t945 ^ 2;
t753 = ((-t856 - t888) * t764 * pkin(3) * t922 + (-(t925 * t978 + t808) * t781 + (-t781 * t856 * t945 + t1103 * t1139) * t1057) * t781 * t1083) * t830;
t756 = (t759 - t1032) * t830 * t781;
t947 = rSges(3,3) + pkin(5);
t1021 = t947 + t1118;
t817 = -t1021 * t1127 + Icges(3,6);
t818 = t1021 * rSges(3,1) * m(3) - Icges(3,5);
t796 = -t817 * t945 + t939 * t818;
t950 = -Icges(3,1) / 0.2e1;
t974 = rSges(3,2) ^ 2;
t975 = rSges(3,1) ^ 2;
t824 = (t975 / 0.2e1 - t974 / 0.2e1) * m(3) + t950 + Icges(3,2) / 0.2e1;
t834 = rSges(3,2) * t1031;
t1059 = t974 + t975;
t844 = -t1059 * m(3) - Icges(3,3);
t866 = -Icges(3,4) + t1047;
t1015 = 0.2e1 * t781 ^ 2 * (t866 * t925 + (t824 * t939 + t834) * t945 + t939 * t1011 + t1063) - t753 * t1119 + t844 * t747 + t796 * t756;
t1138 = t1015 * t790 * t979;
t935 = sin(qJ(3,3));
t941 = cos(qJ(3,3));
t836 = t941 * rSges(3,1) - t935 * rSges(3,2);
t1121 = m(3) * t836;
t886 = t941 * pkin(3);
t861 = t886 + t976;
t1041 = t861 * t879;
t965 = 0.4e1 * qJ(3,3);
t1069 = t978 * cos(t965);
t967 = 0.2e1 * qJ(3,3);
t913 = cos(t967);
t1094 = t816 * t913;
t1096 = t806 * t913;
t1048 = pkin(7) + qJ(3,3);
t1050 = -pkin(7) + qJ(3,3);
t897 = qJ(1,3) + pkin(7);
t875 = cos(t897);
t904 = sin(t967);
t1099 = (t875 * t1129 + sin(t897) * t976 + t936 * t1132 + (sin(qJ(1,3) - t1050) + sin(qJ(1,3) + t1048)) * pkin(3)) / (t935 * t976 + pkin(3) * t904 + (sin(t1048) + sin(t1050)) * pkin(1));
t862 = t886 + pkin(2);
t1091 = t862 * t942;
t786 = (-t1079 + t1091) * t928 - t825 - t862 * t1080;
t828 = 0.1e1 / (t879 + t862);
t932 = legFrame(3,2);
t880 = sin(t932);
t883 = cos(t932);
t920 = 0.1e1 / t935;
t765 = (t952 * t1099 + (t880 * t953 - t883 * t954) * t920 * t828 * t786) * t979;
t1105 = t765 * t935;
t1034 = pkin(3) * t1105;
t850 = t932 + t897;
t851 = -t932 + t897;
t810 = -sin(t850) - sin(t851);
t813 = cos(t851) - cos(t850);
t779 = (t810 * t1143 + t813 * t1144 - t875 * t952) * t828;
t757 = t1139 * t779 - t1034;
t902 = sin(t965);
t966 = 0.3e1 * qJ(3,3);
t903 = sin(t966);
t912 = cos(t966);
t991 = pkin(3) * (-t941 * pkin(2) + t856 * t912);
t746 = ((t912 * t1044 + (t1139 * t861 + t1012 - t1096) * t1057) * t765 + (t903 * t1043 + t904 * t1053 - t977 * t902 + t984 * t935) * t779) / (t1094 + t1069 / 0.2e1 - 0.2e1 * t1041 + 0.2e1 * t991 + t1007) * t779 * t1122 + (t757 * t1131 + (t904 * t1130 - t978 * t902 + (-t856 * t903 - t935 * t879) * t1056) * t765 + (-0.2e1 * t1096 + (-t912 + t941) * t1023 + t995) * t779) / (t1009 - 0.4e1 * t1041 + t1069 + 0.2e1 * t1094 + 0.4e1 * t991) * t765;
t1085 = t920 * t941;
t762 = t765 ^ 2;
t923 = t941 ^ 2;
t751 = ((-t856 - t886) * t762 * pkin(3) * t920 + (-(t923 * t978 + t808) * t779 + (-t779 * t856 * t941 + t1105 * t1139) * t1057) * t779 * t1085) * t828;
t754 = (t757 - t1034) * t828 * t779;
t794 = -t817 * t941 + t935 * t818;
t1016 = 0.2e1 * t779 ^ 2 * (t866 * t923 + (t824 * t935 + t834) * t941 + t935 * t1011 + t1063) - t751 * t1121 + t844 * t746 + t794 * t754;
t1137 = t1016 * t786 * t979;
t838 = t943 * rSges(3,1) - t937 * rSges(3,2);
t1120 = m(3) * t838;
t1010 = t979 * t1026;
t969 = 0.3e1 * qJ(3,2);
t915 = cos(t969);
t1061 = t915 - t943;
t1092 = t829 / 0.2e1;
t852 = t933 + t900;
t853 = -t933 + t900;
t811 = -sin(t852) - sin(t853);
t800 = t811 * t954 * t1092;
t814 = cos(t853) - cos(t852);
t801 = t814 * t953 * t1092;
t1064 = t800 + t801;
t968 = 0.4e1 * qJ(3,2);
t1068 = t978 * cos(t968);
t1088 = t876 * t952;
t766 = (t1026 + (-t1086 + t1087) * t1028) * t979;
t1115 = pkin(3) * t766;
t1104 = t766 * t937;
t1033 = pkin(3) * t1104;
t1025 = t829 * t1088;
t780 = -t1025 + t1064;
t758 = t1139 * t780 - t1033;
t760 = -t1010 / 0.2e1 + (-t1088 + t1140) * t829 + t1064;
t761 = t1010 / 0.2e1 + (-t1088 - t1140) * t829 + t1064;
t775 = t981 * t780;
t858 = 0.2e1 * t893;
t859 = 0.2e1 * t895;
t892 = pkin(7) + t970;
t871 = cos(t892);
t872 = cos(t893);
t894 = -pkin(7) + t970;
t873 = cos(t894);
t874 = cos(t895);
t891 = t964 + qJ(3,2);
t896 = -0.2e1 * pkin(7) + qJ(3,2);
t898 = t969 + pkin(7);
t899 = t969 - pkin(7);
t905 = sin(t968);
t906 = sin(t969);
t916 = cos(t970);
t985 = t775 + (0.3e1 / 0.4e1 * t800 + 0.3e1 / 0.4e1 * t801 - 0.3e1 / 0.4e1 * t1025) * t978 + 0.3e1 * t780 * (t980 + t926 / 0.3e1);
t739 = (t1042 * t1115 - 0.4e1 * ((t985 + t1145) * t870 + (t985 - t1145) * t869) * pkin(1) + 0.4e1 * (t1116 + (pkin(2) * t916 - t879 - t887 / 0.2e1) * t955) * t1115 + 0.4e1 * ((cos(t896) - cos(t891)) * t955 - (sin(t896) + sin(t891)) * pkin(2)) * t775 + ((t760 * sin(t859) + t761 * sin(t858)) * t963 + 0.4e1 * (-t760 * t871 + t761 * t873) * t1112 - 0.12e2 * ((-t1010 / 0.6e1 + (-t1088 + t1142) * t829 + t1064) * sin(t894) + (t1010 / 0.6e1 + (-t1088 - t1142) * t829 + t1064) * sin(t892)) * pkin(1) * pkin(2)) * pkin(3) + (t766 * t915 * t1129 - 0.3e1 * ((-t1010 / 0.3e1 + (-t1088 + t1141) * t829 + t1064) * sin(t899) + (t1010 / 0.3e1 + (-t1088 - t1141) * t829 + t1064) * sin(t898)) * pkin(1)) * t978 + (-0.8e1 * (t978 / 0.4e1 + 0.3e1 / 0.2e1 * t980 + t1062) * t1114 - t977 * t905 - 0.16e2 * t1008 * t1113 + (t906 * t1128 + 0.8e1 * (-t872 + t874) * t1112) * pkin(2)) * t780) / (t959 * t916 + t1068 / 0.2e1 + (cos(t859) / 0.2e1 + cos(t858) / 0.2e1 + t916 - 0.1e1) * t981 + t1061 * pkin(2) * t1057 + ((t871 + t873 - 0.2e1 * t928) * t976 + (cos(t899) + cos(t898) - t872 - t874) * pkin(3)) * pkin(1) + t1024) * t780 * t1122 + (t758 * t1131 + (t907 * t1130 - t978 * t905 + (-t856 * t906 - t937 * t879) * t1056) * t766 + (-t1061 * t1023 - 0.2e1 * t806 * t916 + t995) * t780) / (0.2e1 * t816 * t916 + t1068 - 0.4e1 * (t887 + t976) * t879 + (-t943 * pkin(2) + t856 * t915) * t1056 + t1009) * t766;
t1084 = t921 * t943;
t763 = t766 ^ 2;
t924 = t943 ^ 2;
t752 = ((-t856 - t887) * t763 * pkin(3) * t921 + (-(t924 * t978 + t808) * t780 + (-t780 * t856 * t943 + t1104 * t1139) * t1057) * t780 * t1084) * t829;
t755 = (t758 - t1033) * t829 * t780;
t795 = -t817 * t943 + t937 * t818;
t1018 = 0.2e1 * t780 ^ 2 * (t866 * t924 + (t824 * t937 + t834) * t943 + t937 * t1011 + t1063) - t752 * t1120 + t844 * t739 + t795 * t755;
t1136 = t1018 * t1101;
t1126 = -Icges(3,5) / 0.4e1;
t1125 = Icges(3,6) / 0.4e1;
t843 = (-t974 + t975) * m(3) + Icges(3,2) - Icges(3,1);
t1124 = -t843 / 0.2e1;
t1123 = t843 / 0.2e1;
t1082 = t935 * t880;
t1081 = t935 * t883;
t1078 = t937 * t881;
t1077 = t937 * t884;
t1074 = t939 * t882;
t1073 = t939 * t885;
t1058 = 0.2e1 * m(3);
t1037 = pkin(3) * (-t942 * t928 + t1080) * t923;
t1036 = pkin(3) * (-t944 * t928 + t1076) * t924;
t1035 = pkin(3) * (-t946 * t928 + t1072) * t925;
t1002 = t1118 / 0.4e1 + t947 / 0.4e1;
t983 = -(0.2e1 * t947 ^ 2 + t1059 + t959 + t962) * m(3) / 0.2e1 + 0.2e1 * (rSges(2,2) * m(2) - t947 * m(3)) * t1118 - (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - Icges(3,2) / 0.2e1 + t950 + (-0.2e1 * rSges(2,1) * t879 - t981) * m(2) - Icges(1,3) - Icges(2,3);
t1006 = t879 / 0.2e1 + pkin(2) / 0.2e1;
t997 = t1006 * t780;
t1017 = 0.4e1 * ((t937 * t1125 + t943 * t1126) * t766 + (t937 * t943 * t1123 + (t924 - 0.1e1 / 0.2e1) * t866) * t780 + ((-t1002 * t1104 + t943 * t997) * rSges(3,2) + (t1002 * t943 * t766 + t937 * t997) * rSges(3,1)) * m(3)) * t766 - (t916 * t1124 + t866 * t907 + (-t838 * t879 + (-t879 - t838) * pkin(2)) * t1058 + t983) * t755 - t795 * t739;
t998 = t1006 * t779;
t1014 = 0.4e1 * ((t935 * t1125 + t941 * t1126) * t765 + (t935 * t941 * t1123 + (t923 - 0.1e1 / 0.2e1) * t866) * t779 + ((-t1002 * t1105 + t941 * t998) * rSges(3,2) + (t1002 * t941 * t765 + t935 * t998) * rSges(3,1)) * m(3)) * t765 - (t913 * t1124 + t866 * t904 + (-t836 * t879 + (-t879 - t836) * pkin(2)) * t1058 + t983) * t754 - t794 * t746;
t996 = t1006 * t781;
t1013 = 0.4e1 * ((t939 * t1125 + t945 * t1126) * t767 + (t939 * t945 * t1123 + (t925 - 0.1e1 / 0.2e1) * t866) * t781 + ((-t1002 * t1103 + t945 * t996) * rSges(3,2) + (t1002 * t945 * t767 + t939 * t996) * rSges(3,1)) * m(3)) * t767 - (t919 * t1124 + t866 * t910 + (-t840 * t879 + (-t879 - t840) * pkin(2)) * t1058 + t983) * t756 - t796 * t747;
t957 = -m(2) - m(3);
t1005 = m(3) * t763 * (t937 * rSges(3,1) + t943 * rSges(3,2)) + t739 * t1120 - t957 * t752;
t1004 = m(3) * t762 * (t935 * rSges(3,1) + t941 * rSges(3,2)) + t746 * t1121 - t957 * t751;
t1003 = m(3) * t764 * (t939 * rSges(3,1) + t945 * rSges(3,2)) + t747 * t1119 - t957 * t753;
t1001 = -t1017 / 0.2e1;
t1000 = -t1014 / 0.2e1;
t999 = -t1013 / 0.2e1;
t1 = [(t812 * t999 + (-t885 * t1138 - t1003 * (-t885 * t1035 + (pkin(3) * t1074 - t1146 * t885) * t945 + t856 * t1074)) * t922) * t830 + (t811 * t1001 + (-t884 * t1136 - t1005 * (-t884 * t1036 + (pkin(3) * t1078 - t1147 * t884) * t943 + t856 * t1078)) * t921) * t829 + (t810 * t1000 + (-t883 * t1137 - t1004 * (-t883 * t1037 + (pkin(3) * t1082 - t1148 * t883) * t941 + t856 * t1082)) * t920) * t828; (t815 * t999 + (t882 * t1138 - t1003 * (t882 * t1035 + (pkin(3) * t1073 + t1146 * t882) * t945 + t856 * t1073)) * t922) * t830 + (t814 * t1001 + (t881 * t1136 - t1005 * (t881 * t1036 + (pkin(3) * t1077 + t1147 * t881) * t943 + t856 * t1077)) * t921) * t829 + (t813 * t1000 + (t880 * t1137 - t1004 * (t880 * t1037 + (pkin(3) * t1081 + t1148 * t880) * t941 + t856 * t1081)) * t920) * t828; (t1013 * t877 + t1003 * ((t865 * t940 + t946 * t955) * t928 - t845 * t940 + t927 * t1089) * t1083) * t830 + (t1017 * t876 + t1005 * ((t863 * t938 + t944 * t955) * t928 - t845 * t938 + t927 * t1090) * t1084) * t829 + (t1014 * t875 + t1004 * ((t862 * t936 + t942 * t955) * t928 - t845 * t936 + t927 * t1091) * t1085) * t828 + (t1015 * t1097 + t1016 * t1099 + t1018 * t1098) * t979;];
taucX  = t1;
