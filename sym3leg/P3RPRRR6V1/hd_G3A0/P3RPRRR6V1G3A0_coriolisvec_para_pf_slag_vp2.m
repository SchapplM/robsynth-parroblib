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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:41:45
% EndTime: 2020-08-06 18:41:53
% DurationCPUTime: 7.29s
% Computational Cost: add. (43605->450), mult. (48396->644), div. (7563->16), fcn. (33165->88), ass. (0->333)
t951 = sin(pkin(7));
t1150 = t951 * pkin(1);
t1175 = pkin(5) + t1150;
t1184 = t1175 * mrSges(3,2);
t952 = cos(pkin(7));
t896 = t952 * pkin(1);
t868 = t896 + pkin(2);
t1183 = mrSges(3,1) * t868;
t965 = sin(qJ(1,1));
t978 = -pkin(6) - pkin(5);
t1108 = t965 * t978;
t1109 = t965 * t951;
t1107 = t978 * t951;
t853 = -pkin(1) + t1107;
t971 = cos(qJ(1,1));
t836 = t853 * t971;
t1182 = -pkin(2) * t1109 + (pkin(2) * t971 - t1108) * t952 - t836;
t963 = sin(qJ(1,2));
t1112 = t963 * t978;
t1113 = t963 * t951;
t969 = cos(qJ(1,2));
t835 = t853 * t969;
t1181 = -pkin(2) * t1113 + (pkin(2) * t969 - t1112) * t952 - t835;
t961 = sin(qJ(1,3));
t1116 = t961 * t978;
t1117 = t961 * t951;
t967 = cos(qJ(1,3));
t834 = t853 * t967;
t1180 = -pkin(2) * t1117 + (pkin(2) * t967 - t1116) * t952 - t834;
t1159 = -2 * Ifges(3,4);
t1176 = mrSges(3,2) * t868;
t970 = cos(qJ(3,1));
t949 = t970 ^ 2;
t956 = Ifges(3,1) - Ifges(3,2);
t964 = sin(qJ(3,1));
t1179 = t949 * t1159 + Ifges(3,4) + (-t956 * t964 + t1176) * t970;
t968 = cos(qJ(3,2));
t948 = t968 ^ 2;
t962 = sin(qJ(3,2));
t1178 = t948 * t1159 + Ifges(3,4) + (-t956 * t962 + t1176) * t968;
t966 = cos(qJ(3,3));
t947 = t966 ^ 2;
t960 = sin(qJ(3,3));
t1177 = t947 * t1159 + Ifges(3,4) + (-t956 * t960 + t1176) * t966;
t997 = 0.2e1 * pkin(2);
t1174 = 0.4e1 * pkin(2);
t976 = xDP(2);
t1172 = t976 / 0.2e1;
t977 = xDP(1);
t1171 = t977 / 0.2e1;
t1001 = pkin(2) ^ 2;
t1002 = pkin(1) ^ 2;
t980 = m(2) + m(3);
t1006 = -(pkin(5) ^ 2 + t1001) * m(3) - t980 * t1002 - Ifges(3,1) - Ifges(1,3) - Ifges(2,3) - 0.2e1 * (m(3) * pkin(5) - mrSges(2,2) + mrSges(3,3)) * t1150 - 0.2e1 * mrSges(3,3) * pkin(5);
t1077 = mrSges(3,1) * t1150;
t898 = mrSges(3,1) * pkin(5) - Ifges(3,5);
t1018 = t1077 / 0.2e1 + t898 / 0.2e1;
t1019 = Ifges(3,6) / 0.2e1 - t1184 / 0.2e1;
t1062 = -m(3) * pkin(2) - mrSges(2,1);
t1149 = t962 * pkin(2);
t1086 = 0.2e1 * t1149;
t1088 = 0.2e1 * t896;
t1000 = 0.1e1 / pkin(3);
t909 = t968 * pkin(3);
t875 = t909 + pkin(2);
t1128 = t875 * t969;
t799 = (-t1112 + t1128) * t952 - t835 - t875 * t1113;
t958 = legFrame(2,2);
t906 = cos(t958);
t945 = 0.1e1 / t962;
t1060 = t799 * t906 * t945;
t1021 = t1000 * t977 * t1060;
t838 = 0.1e1 / (t896 + t875);
t1007 = t838 * t1021;
t1132 = t838 * t945;
t903 = sin(t958);
t1031 = t799 * t903 * t1132;
t794 = t1000 * t976 * t1031;
t993 = 0.2e1 * qJ(3,2);
t930 = sin(t993);
t1151 = pkin(3) * t930;
t1157 = 0.2e1 * t978;
t1162 = 0.2e1 * pkin(1);
t916 = pkin(7) + qJ(3,2);
t881 = sin(t916);
t918 = -pkin(7) + qJ(3,2);
t883 = sin(t918);
t923 = qJ(1,2) + pkin(7);
t893 = cos(t923);
t1139 = (t893 * t1157 + sin(t923) * t997 + t963 * t1162 + (sin(qJ(1,2) - t918) + sin(qJ(1,2) + t916)) * pkin(3)) / (t1086 + t1151 + (t881 + t883) * pkin(1));
t975 = xDP(3);
t795 = t1000 * t975 * t1139;
t1004 = t794 / 0.4e1 + t795 / 0.4e1 - t1007 / 0.4e1;
t887 = cos(t916);
t889 = cos(t918);
t992 = 0.3e1 * qJ(3,2);
t921 = t992 + pkin(7);
t922 = t992 - pkin(7);
t1020 = (-t887 - t889 + cos(t921) + cos(t922)) * pkin(3);
t915 = pkin(7) + t993;
t886 = cos(t915);
t917 = -pkin(7) + t993;
t888 = cos(t917);
t1056 = t886 + t888 - 0.2e1 * t952;
t1023 = 0.2e1 * t1056;
t938 = cos(t992);
t1049 = pkin(3) * (t938 - t968);
t1026 = 0.2e1 * t1049;
t950 = t978 ^ 2;
t1097 = t1002 + t950 / 0.2e1;
t999 = pkin(3) ^ 2;
t1029 = 0.3e1 / 0.8e1 * t999 + t1001 / 0.2e1 + t1097;
t987 = 0.2e1 * pkin(7);
t913 = cos(t987);
t879 = t1002 * t913;
t986 = -0.2e1 * t1002;
t1030 = -0.2e1 * t879 - 0.4e1 * t1001 + t986 - t999;
t1050 = -t794 / 0.6e1 - t795 / 0.6e1;
t1051 = -t794 / 0.3e1 - t795 / 0.3e1;
t1164 = -t879 - t999 / 0.2e1;
t1055 = -0.2e1 * t1001 + t1164;
t1126 = t893 * t975;
t1058 = t838 * t1126;
t1131 = t838 / 0.2e1;
t864 = t958 + t923;
t865 = -t958 + t923;
t825 = -sin(t864) - sin(t865);
t811 = t825 * t977 * t1131;
t828 = cos(t865) - cos(t864);
t812 = t828 * t976 * t1131;
t1057 = 0.3e1 / 0.4e1 * t811 + 0.3e1 / 0.4e1 * t812 - 0.3e1 / 0.4e1 * t1058;
t912 = sin(t987);
t1098 = t1002 * t912;
t1064 = 0.2e1 * t1098;
t982 = 0.2e1 * t1001;
t1095 = t982 + t1002;
t1100 = t811 + t812;
t1101 = t794 + t795;
t792 = -t1058 + t1100;
t787 = t1002 * t792;
t1103 = t787 + 0.3e1 * (t1001 + t950 / 0.3e1) * t792;
t991 = 0.4e1 * qJ(3,2);
t1105 = t999 * cos(t991);
t1141 = t792 * t978;
t1148 = t1000 / 0.2e1;
t778 = -t1007 + t1101;
t1152 = pkin(3) * t778;
t1153 = pkin(2) * sin(t992);
t861 = pkin(6) + t1175;
t1154 = pkin(2) * t861;
t1156 = -0.6e1 * t999;
t1160 = -0.4e1 * pkin(3);
t1165 = 0.4e1 * t978;
t1044 = -t978 + t1150;
t1143 = t778 * t962;
t1066 = pkin(3) * t1143;
t769 = t1044 * t792 - t1066;
t1052 = -t794 / 0.2e1 - t795 / 0.2e1;
t770 = (-t1126 + t1021 / 0.2e1) * t838 + t1052 + t1100;
t771 = (-t1126 - t1021 / 0.2e1) * t838 - t1052 + t1100;
t1102 = 0.2e1 * t1101;
t772 = (0.2e1 * t1021 - t1126) * t838 + t1100 - t1102;
t773 = (-0.2e1 * t1021 - t1126) * t838 + t1100 + t1102;
t870 = 0.2e1 * t916;
t857 = sin(t870);
t871 = 0.2e1 * t918;
t858 = sin(t871);
t859 = cos(t870);
t860 = cos(t871);
t880 = sin(t915);
t882 = sin(t917);
t884 = sin(t921);
t885 = sin(t922);
t914 = t987 + qJ(3,2);
t919 = -0.2e1 * pkin(7) + qJ(3,2);
t928 = sin(t991);
t939 = cos(t993);
t985 = 0.2e1 * t1002;
t998 = pkin(3) * t999;
t749 = (t999 * t778 * t938 * t1157 + t1064 * t1152 + 0.4e1 * (t1154 + (pkin(2) * t939 - t896 - t909 / 0.2e1) * t978) * t1152 + (t770 * t858 + t771 * t857) * pkin(3) * t986 + (-0.8e1 * (t999 / 0.4e1 + 0.3e1 / 0.2e1 * t1001 + t1097) * t1151 + t1153 * t1156 - t998 * t928 - 0.16e2 * t1029 * t1149) * t792 + ((cos(t919) - cos(t914)) * t1165 - 0.4e1 * (sin(t919) + sin(t914)) * pkin(2)) * t787 + (((-t770 * t886 + t771 * t888) * t1165 - 0.12e2 * (((-t1126 + t1021 / 0.6e1) * t838 + t1050 + t1100) * t882 + ((-t1126 - t1021 / 0.6e1) * t838 - t1050 + t1100) * t880) * pkin(2)) * pkin(3) - 0.4e1 * ((t1004 + t1057) * t999 + t1103) * t883 - 0.4e1 * ((-t1004 + t1057) * t999 + t1103) * t881 - 0.3e1 * (((-t1126 + t1021 / 0.3e1) * t838 + t1051 + t1100) * t885 + ((-t1126 - t1021 / 0.3e1) * t838 - t1051 + t1100) * t884) * t999 + 0.8e1 * (-t887 + t889) * pkin(2) * t1141) * pkin(1)) / (t982 * t939 + t1105 / 0.2e1 + (t860 / 0.2e1 + t859 / 0.2e1 + t939 - 0.1e1) * t1002 + pkin(2) * t1026 + (pkin(2) * t1023 + t1020) * pkin(1) + t1055) * t792 * t1148 + (t769 * t1174 + ((t1100 - t1101) * t858 - (t1100 + t1101) * t857 + ((t1021 - t1126) * t858 - (-t1021 - t1126) * t857) * t838) * t1002 + (-0.2e1 * (t999 + t1095) * t930 + t1153 * t1160 - t999 * t928) * t778 + (t1064 + (t1174 * t939 + t1026) * t978) * t792 + ((t772 * t882 - t773 * t880) * t997 + t1023 * t1141 + ((-t883 - t884) * t773 + (t881 + t885) * t772) * pkin(3)) * pkin(1)) / ((t985 + 0.4e1 * t1001) * t939 + t1105 + (t860 + t859) * t1002 + t1049 * t1174 + (t1056 * t1174 + 0.2e1 * t1020) * pkin(1) + t1030) * t778;
t766 = (t769 - t1066) * t838 * t792;
t844 = t1077 + t898;
t845 = Ifges(3,6) - t1184;
t814 = t962 * t844 - t845 * t968;
t900 = t962 * mrSges(3,2);
t1039 = 0.2e1 * ((t1018 * t968 + t1019 * t962) * t778 + (t1183 * t962 + t1178) * t792) * t778 - (t956 * t948 - 0.2e1 * (Ifges(3,4) * t962 + t1183) * t968 + (t1062 + t900) * t1088 + mrSges(3,2) * t1086 + t1006) * t766 - t814 * t749;
t1170 = -t1039 * t838 / 0.2e1;
t1085 = t964 * t997;
t1072 = pkin(1) * t1107;
t1096 = t1002 + t950;
t981 = 0.3e1 * t1001;
t1005 = -0.8e1 * (0.3e1 / 0.4e1 * t999 + t981 + t1096) * t896 - 0.8e1 * (t1029 - t1072) * t997 + 0.8e1 * (-pkin(2) * t913 + t912 * t978) * t1002;
t995 = 0.3e1 * qJ(3,1);
t941 = cos(t995);
t1011 = pkin(3) * (-t970 * pkin(2) + t868 * t941);
t1071 = t978 * t896;
t1013 = t1064 - 0.4e1 * t1071;
t1028 = -t1002 + t1055;
t1032 = -0.2e1 * t1071 + t1098;
t1093 = 0.2e1 * pkin(3);
t1048 = t861 * t1093;
t910 = t970 * pkin(3);
t876 = t910 + t997;
t1073 = t876 * t896;
t1075 = t868 * t1156;
t1076 = -0.2e1 * t861 * t999;
t1079 = pkin(2) * t896;
t855 = -0.2e1 * t1072;
t1089 = (0.6e1 * t1079 + t855 + t985 + t981 + t950 - t1164) * t1160;
t1092 = 0.4e1 * pkin(3);
t994 = 0.4e1 * qJ(3,1);
t1104 = t999 * cos(t994);
t830 = 0.4e1 * t1079 + t879 + t1095;
t996 = 0.2e1 * qJ(3,1);
t942 = cos(t996);
t1134 = t830 * t942;
t820 = t1032 + 0.2e1 * t1154;
t1136 = t820 * t942;
t1158 = -0.2e1 * t999 - 0.2e1 * t830;
t877 = t910 + pkin(2);
t1127 = t877 * t971;
t839 = 0.1e1 / (t896 + t877);
t946 = 0.1e1 / t964;
t1130 = t839 * t946;
t1059 = ((-t1108 + t1127) * t952 - t836 - t877 * t1109) * t1130;
t1082 = pkin(7) + qJ(3,1);
t1084 = -pkin(7) + qJ(3,1);
t924 = qJ(1,1) + pkin(7);
t894 = cos(t924);
t933 = sin(t996);
t1138 = (t894 * t1157 + sin(t924) * t997 + t965 * t1162 + (sin(qJ(1,1) - t1084) + sin(qJ(1,1) + t1082)) * pkin(3)) / (t1085 + pkin(3) * t933 + (sin(t1082) + sin(t1084)) * pkin(1));
t959 = legFrame(1,2);
t904 = sin(t959);
t907 = cos(t959);
t779 = (t975 * t1138 + (t904 * t976 - t907 * t977) * t1059) * t1000;
t1142 = t779 * t964;
t1065 = pkin(3) * t1142;
t866 = t959 + t924;
t867 = -t959 + t924;
t826 = -sin(t866) - sin(t867);
t829 = cos(t867) - cos(t866);
t793 = (t826 * t1171 + t829 * t1172 - t894 * t975) * t839;
t767 = t1044 * t793 - t1065;
t931 = sin(t994);
t932 = sin(t995);
t757 = ((t941 * t1076 + (t876 * t861 + t1032 - t1136) * t1093) * t779 + (t1005 * t964 + t932 * t1075 + t933 * t1089 - t998 * t931) * t793) / (t1134 + t1104 / 0.2e1 - 0.2e1 * t1073 + 0.2e1 * t1011 + t1028) * t793 * t1148 + (t767 * t1174 + (t933 * t1158 - t999 * t931 + (-t868 * t932 - t964 * t896) * t1092) * t779 + (-0.2e1 * t1136 + (-t941 + t970) * t1048 + t1013) * t793) / (0.4e1 * t1011 + t1030 - 0.4e1 * t1073 + t1104 + 0.2e1 * t1134) * t779;
t764 = (t767 - t1065) * t839 * t793;
t815 = t964 * t844 - t845 * t970;
t901 = t964 * mrSges(3,2);
t1033 = 0.2e1 * ((t1018 * t970 + t1019 * t964) * t779 + (t1183 * t964 + t1179) * t793) * t779 - (t956 * t949 - 0.2e1 * (Ifges(3,4) * t964 + t1183) * t970 + (t1062 + t901) * t1088 + mrSges(3,2) * t1085 + t1006) * t764 - t815 * t757;
t1169 = -t839 * t1033 / 0.2e1;
t1087 = t960 * t997;
t989 = 0.3e1 * qJ(3,3);
t935 = cos(t989);
t1012 = pkin(3) * (-t966 * pkin(2) + t868 * t935);
t908 = t966 * pkin(3);
t873 = t908 + t997;
t1074 = t873 * t896;
t988 = 0.4e1 * qJ(3,3);
t1106 = t999 * cos(t988);
t990 = 0.2e1 * qJ(3,3);
t936 = cos(t990);
t1135 = t830 * t936;
t1137 = t820 * t936;
t874 = t908 + pkin(2);
t1129 = t874 * t967;
t837 = 0.1e1 / (t896 + t874);
t944 = 0.1e1 / t960;
t1133 = t837 * t944;
t1061 = ((-t1116 + t1129) * t952 - t834 - t874 * t1117) * t1133;
t1081 = pkin(7) + qJ(3,3);
t1083 = -pkin(7) + qJ(3,3);
t920 = qJ(1,3) + pkin(7);
t890 = cos(t920);
t927 = sin(t990);
t1140 = (t890 * t1157 + sin(t920) * t997 + t961 * t1162 + (sin(qJ(1,3) - t1083) + sin(qJ(1,3) + t1081)) * pkin(3)) / (t1087 + pkin(3) * t927 + (sin(t1081) + sin(t1083)) * pkin(1));
t957 = legFrame(3,2);
t902 = sin(t957);
t905 = cos(t957);
t777 = (t975 * t1140 + (t902 * t976 - t905 * t977) * t1061) * t1000;
t1144 = t777 * t960;
t1067 = pkin(3) * t1144;
t862 = t957 + t920;
t863 = -t957 + t920;
t824 = -sin(t862) - sin(t863);
t827 = cos(t863) - cos(t862);
t791 = (t824 * t1171 + t827 * t1172 - t890 * t975) * t837;
t768 = t1044 * t791 - t1067;
t925 = sin(t988);
t926 = sin(t989);
t756 = ((t935 * t1076 + (t873 * t861 + t1032 - t1137) * t1093) * t777 + (t1005 * t960 + t926 * t1075 + t927 * t1089 - t998 * t925) * t791) / (t1135 + t1106 / 0.2e1 - 0.2e1 * t1074 + 0.2e1 * t1012 + t1028) * t791 * t1148 + (t768 * t1174 + (t927 * t1158 - t999 * t925 + (-t868 * t926 - t960 * t896) * t1092) * t777 + (-0.2e1 * t1137 + (-t935 + t966) * t1048 + t1013) * t791) / (0.4e1 * t1012 + t1030 - 0.4e1 * t1074 + t1106 + 0.2e1 * t1135) * t777;
t765 = (t768 - t1067) * t837 * t791;
t813 = t960 * t844 - t845 * t966;
t899 = t960 * mrSges(3,2);
t1034 = 0.2e1 * ((t1018 * t966 + t1019 * t960) * t777 + (t1183 * t960 + t1177) * t791) * t777 - (t956 * t947 - 0.2e1 * (Ifges(3,4) * t960 + t1183) * t966 + (t1062 + t899) * t1088 + mrSges(3,2) * t1087 + t1006) * t765 - t813 * t756;
t1168 = -t837 * t1034 / 0.2e1;
t1145 = mrSges(3,1) * t964;
t1123 = t946 * t970;
t776 = t779 ^ 2;
t822 = t1001 + t855 + 0.2e1 * t1079 + t1096;
t763 = (pkin(3) * (-t868 - t910) * t776 * t946 + (-(t949 * t999 + t822) * t793 + (-t793 * t868 * t970 + t861 * t1142) * t1093) * t793 * t1123) * t839;
t846 = -t970 * mrSges(3,1) + t901;
t1037 = t793 ^ 2 * (t868 * t1145 + t1179) - Ifges(3,3) * t757 + t846 * t763 + t815 * t764;
t1167 = t1037 * t1059;
t1147 = mrSges(3,1) * t960;
t1125 = t944 * t966;
t774 = t777 ^ 2;
t761 = (pkin(3) * (-t868 - t908) * t774 * t944 + (-(t947 * t999 + t822) * t791 + (-t791 * t868 * t966 + t861 * t1144) * t1093) * t791 * t1125) * t837;
t847 = -t966 * mrSges(3,1) + t899;
t1038 = t791 ^ 2 * (t868 * t1147 + t1177) - Ifges(3,3) * t756 + t847 * t761 + t813 * t765;
t1166 = t1038 * t1061;
t775 = t778 ^ 2;
t1146 = mrSges(3,1) * t962;
t1124 = t945 * t968;
t1119 = t960 * t902;
t1118 = t960 * t905;
t1115 = t962 * t903;
t1114 = t962 * t906;
t1111 = t964 * t904;
t1110 = t964 * t907;
t1070 = pkin(3) * (-t967 * t952 + t1117) * t947;
t1069 = pkin(3) * (-t969 * t952 + t1113) * t948;
t1068 = pkin(3) * (-t971 * t952 + t1109) * t949;
t762 = (pkin(3) * (-t868 - t909) * t775 * t945 + (-(t948 * t999 + t822) * t792 + (-t792 * t868 * t968 + t861 * t1143) * t1093) * t792 * t1124) * t838;
t848 = -t968 * mrSges(3,1) + t900;
t1041 = t792 ^ 2 * (t868 * t1146 + t1178) - Ifges(3,3) * t749 + t848 * t762 + t814 * t766;
t1040 = -t848 * t749 + t980 * t762 + t775 * (mrSges(3,2) * t968 + t1146);
t1036 = -t847 * t756 + t980 * t761 + t774 * (mrSges(3,2) * t966 + t1147);
t1035 = -t846 * t757 + t980 * t763 + t776 * (mrSges(3,2) * t970 + t1145);
t1010 = t1036 * t1133;
t1009 = t1040 * t1132;
t1008 = t1035 * t1130;
t1 = [-(-t907 * t1068 + (pkin(3) * t1111 + t1182 * t907) * t970 + t868 * t1111) * t1008 - (-t906 * t1069 + (pkin(3) * t1115 + t1181 * t906) * t968 + t868 * t1115) * t1009 - (-t905 * t1070 + (pkin(3) * t1119 + t1180 * t905) * t966 + t868 * t1119) * t1010 + t825 * t1170 + t826 * t1169 + t824 * t1168 + (-t1041 * t838 * t1060 - t905 * t1166 - t907 * t1167) * t1000; -(t904 * t1068 + (pkin(3) * t1110 - t1182 * t904) * t970 + t868 * t1110) * t1008 - (t903 * t1069 + (pkin(3) * t1114 - t1181 * t903) * t968 + t868 * t1114) * t1009 - (t902 * t1070 + (pkin(3) * t1118 - t1180 * t902) * t966 + t868 * t1118) * t1010 + t828 * t1170 + t829 * t1169 + t827 * t1168 + (t1041 * t1031 + t902 * t1166 + t904 * t1167) * t1000; (t1033 * t894 + t1035 * ((t877 * t965 + t971 * t978) * t952 - t853 * t965 + t951 * t1127) * t1123) * t839 + (t1039 * t893 + t1040 * ((t875 * t963 + t969 * t978) * t952 - t853 * t963 + t951 * t1128) * t1124) * t838 + (t1034 * t890 + t1036 * ((t874 * t961 + t967 * t978) * t952 - t853 * t961 + t951 * t1129) * t1125) * t837 + (t1037 * t1138 + t1038 * t1140 + t1041 * t1139) * t1000;];
taucX  = t1;
