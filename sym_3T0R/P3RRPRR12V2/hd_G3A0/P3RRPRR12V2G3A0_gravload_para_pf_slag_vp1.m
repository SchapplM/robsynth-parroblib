% Calculate Gravitation load for parallel robot
% P3RRPRR12V2G3A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V2G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:27:08
% EndTime: 2020-08-06 19:27:09
% DurationCPUTime: 0.93s
% Computational Cost: add. (870->181), mult. (1173->319), div. (45->6), fcn. (648->18), ass. (0->129)
t1113 = (pkin(2) + pkin(3));
t1177 = -2 * t1113;
t1096 = legFrame(3,2);
t1078 = sin(t1096);
t1176 = t1078 * qJ(3,3);
t1097 = legFrame(2,2);
t1079 = sin(t1097);
t1175 = t1079 * qJ(3,2);
t1098 = legFrame(1,2);
t1080 = sin(t1098);
t1174 = t1080 * qJ(3,1);
t1081 = cos(t1096);
t1173 = t1081 * qJ(3,3);
t1082 = cos(t1097);
t1172 = t1082 * qJ(3,2);
t1083 = cos(t1098);
t1171 = t1083 * qJ(3,1);
t1099 = sin(qJ(2,3));
t1075 = t1099 * qJ(3,3);
t1101 = sin(qJ(2,2));
t1076 = t1101 * qJ(3,2);
t1103 = sin(qJ(2,1));
t1077 = t1103 * qJ(3,1);
t1056 = t1078 * g(1) + t1081 * g(2);
t1111 = m(2) * rSges(2,2);
t1065 = (-qJ(3,3) - rSges(3,3)) * m(3) + t1111;
t1068 = (pkin(2) + rSges(3,1)) * m(3) + m(2) * rSges(2,1);
t1105 = cos(qJ(2,3));
t1059 = t1081 * g(1) - t1078 * g(2);
t1100 = sin(qJ(1,3));
t1106 = cos(qJ(1,3));
t1122 = -g(3) * t1100 + t1059 * t1106;
t1035 = (-t1068 * t1056 + t1122 * t1065) * t1105 + (t1065 * t1056 + t1122 * t1068) * t1099;
t1114 = 0.1e1 / qJ(3,3);
t1170 = t1035 * t1114;
t1057 = t1079 * g(1) + t1082 * g(2);
t1066 = (-qJ(3,2) - rSges(3,3)) * m(3) + t1111;
t1107 = cos(qJ(2,2));
t1060 = t1082 * g(1) - t1079 * g(2);
t1102 = sin(qJ(1,2));
t1108 = cos(qJ(1,2));
t1121 = -g(3) * t1102 + t1060 * t1108;
t1036 = (-t1068 * t1057 + t1121 * t1066) * t1107 + (t1066 * t1057 + t1121 * t1068) * t1101;
t1115 = 0.1e1 / qJ(3,2);
t1169 = t1036 * t1115;
t1058 = t1080 * g(1) + t1083 * g(2);
t1067 = (-qJ(3,1) - rSges(3,3)) * m(3) + t1111;
t1109 = cos(qJ(2,1));
t1061 = t1083 * g(1) - t1080 * g(2);
t1104 = sin(qJ(1,1));
t1110 = cos(qJ(1,1));
t1120 = -g(3) * t1104 + t1061 * t1110;
t1037 = (-t1068 * t1058 + t1120 * t1067) * t1109 + (t1067 * t1058 + t1120 * t1068) * t1103;
t1116 = 0.1e1 / qJ(3,1);
t1168 = t1037 * t1116;
t1132 = t1106 * t1075;
t1112 = pkin(5) - pkin(6);
t1140 = pkin(1) * t1106 + t1100 * t1112;
t1167 = (0.2e1 * t1132 + t1140) * t1113;
t1133 = t1108 * t1076;
t1139 = pkin(1) * t1108 + t1102 * t1112;
t1166 = (0.2e1 * t1133 + t1139) * t1113;
t1134 = t1110 * t1077;
t1138 = pkin(1) * t1110 + t1104 * t1112;
t1165 = (0.2e1 * t1134 + t1138) * t1113;
t1143 = t1113 * t1105;
t1125 = t1075 + pkin(1) + t1143;
t1053 = 0.1e1 / t1125;
t1164 = t1053 * t1114;
t1142 = t1113 * t1107;
t1124 = t1076 + pkin(1) + t1142;
t1054 = 0.1e1 / t1124;
t1163 = t1054 * t1115;
t1141 = t1113 * t1109;
t1123 = t1077 + pkin(1) + t1141;
t1055 = 0.1e1 / t1123;
t1162 = t1055 * t1116;
t1161 = (qJ(3,3) + t1113) * (-qJ(3,3) + t1113);
t1160 = (qJ(3,2) + t1113) * (-qJ(3,2) + t1113);
t1159 = (qJ(3,1) + t1113) * (-qJ(3,1) + t1113);
t1158 = t1099 * t1113;
t1062 = (-rSges(3,2) - pkin(5)) * m(3) + (-rSges(2,3) - pkin(5)) * m(2) + m(1) * rSges(1,2);
t1064 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t1063 = g(3) * t1064;
t1119 = -t1065 * t1099 + t1068 * t1105;
t1032 = t1063 * t1106 + (-t1062 * t1100 + t1119 * t1106) * g(3) + (t1062 * t1106 + (t1064 + t1119) * t1100) * t1059;
t1157 = t1100 * t1032;
t1156 = t1101 * t1113;
t1118 = -t1066 * t1101 + t1068 * t1107;
t1033 = t1063 * t1108 + (-t1062 * t1102 + t1118 * t1108) * g(3) + (t1062 * t1108 + (t1064 + t1118) * t1102) * t1060;
t1155 = t1102 * t1033;
t1154 = t1103 * t1113;
t1117 = -t1067 * t1103 + t1068 * t1109;
t1034 = t1063 * t1110 + (-t1062 * t1104 + t1117 * t1110) * g(3) + (t1062 * t1110 + (t1064 + t1117) * t1104) * t1061;
t1153 = t1104 * t1034;
t1152 = t1106 * t1113;
t1151 = t1108 * t1113;
t1150 = t1110 * t1113;
t1149 = t1112 * t1106;
t1148 = t1112 * t1108;
t1147 = t1112 * t1110;
t1069 = pkin(1) * t1099 + qJ(3,3);
t1146 = t1113 * t1069;
t1070 = pkin(1) * t1101 + qJ(3,2);
t1145 = t1113 * t1070;
t1071 = pkin(1) * t1103 + qJ(3,1);
t1144 = t1113 * t1071;
t1137 = qJ(3,1) * t1177;
t1136 = qJ(3,2) * t1177;
t1135 = qJ(3,3) * t1177;
t1131 = (-t1056 * t1105 + t1122 * t1099) * t1164;
t1130 = (-t1057 * t1107 + t1121 * t1101) * t1163;
t1129 = (-t1058 * t1109 + t1120 * t1103) * t1162;
t1128 = t1106 * t1161;
t1127 = t1108 * t1160;
t1126 = t1110 * t1159;
t1095 = t1109 ^ 2;
t1094 = t1107 ^ 2;
t1093 = t1105 ^ 2;
t1052 = pkin(1) * qJ(3,1) - t1103 * t1159;
t1051 = pkin(1) * qJ(3,2) - t1101 * t1160;
t1050 = pkin(1) * qJ(3,3) - t1099 * t1161;
t1049 = t1134 + t1138;
t1047 = t1133 + t1139;
t1045 = t1132 + t1140;
t1043 = qJ(3,1) * t1110 + t1138 * t1103;
t1042 = qJ(3,2) * t1108 + t1139 * t1101;
t1041 = qJ(3,3) * t1106 + t1140 * t1099;
t1 = [-m(4) * g(1) + (-t1083 * t1153 + ((t1083 * t1150 - t1174) * t1095 + (t1049 * t1083 + t1080 * t1154) * t1109 + t1080 * t1071) * t1168) * t1055 + (-t1082 * t1155 + ((t1082 * t1151 - t1175) * t1094 + (t1047 * t1082 + t1079 * t1156) * t1107 + t1079 * t1070) * t1169) * t1054 + (-t1081 * t1157 + ((t1081 * t1152 - t1176) * t1093 + (t1045 * t1081 + t1078 * t1158) * t1105 + t1078 * t1069) * t1170) * t1053 + (-((t1080 * t1137 + t1083 * t1126) * t1095 + (-t1052 * t1080 + t1083 * t1165) * t1109 + t1043 * t1171 + t1080 * t1144) * t1129 - ((t1079 * t1136 + t1082 * t1127) * t1094 + (-t1051 * t1079 + t1082 * t1166) * t1107 + t1042 * t1172 + t1079 * t1145) * t1130 - ((t1078 * t1135 + t1081 * t1128) * t1093 + (-t1050 * t1078 + t1081 * t1167) * t1105 + t1041 * t1173 + t1078 * t1146) * t1131) * m(3); -m(4) * g(2) + (t1080 * t1153 + ((-t1080 * t1150 - t1171) * t1095 + (-t1049 * t1080 + t1083 * t1154) * t1109 + t1083 * t1071) * t1168) * t1055 + (t1079 * t1155 + ((-t1079 * t1151 - t1172) * t1094 + (-t1047 * t1079 + t1082 * t1156) * t1107 + t1082 * t1070) * t1169) * t1054 + (t1078 * t1157 + ((-t1078 * t1152 - t1173) * t1093 + (-t1045 * t1078 + t1081 * t1158) * t1105 + t1081 * t1069) * t1170) * t1053 + (-((-t1080 * t1126 + t1083 * t1137) * t1095 + (-t1052 * t1083 - t1080 * t1165) * t1109 - t1043 * t1174 + t1083 * t1144) * t1129 - ((-t1079 * t1127 + t1082 * t1136) * t1094 + (-t1051 * t1082 - t1079 * t1166) * t1107 - t1042 * t1175 + t1082 * t1145) * t1130 - ((-t1078 * t1128 + t1081 * t1135) * t1093 + (-t1050 * t1081 - t1078 * t1167) * t1105 - t1041 * t1176 + t1081 * t1146) * t1131) * m(3); -t1110 * t1055 * t1034 - t1109 * (t1123 * t1104 - t1147) * t1037 * t1162 - t1108 * t1054 * t1033 - t1107 * (t1124 * t1102 - t1148) * t1036 * t1163 - t1106 * t1053 * t1032 - t1105 * (t1125 * t1100 - t1149) * t1035 * t1164 - m(4) * g(3) + (-(-t1104 * t1095 * t1159 - ((0.2e1 * t1077 + pkin(1)) * t1104 - t1147) * t1141 - qJ(3,1) * (t1071 * t1104 - t1103 * t1147)) * t1129 - (-t1102 * t1094 * t1160 - ((0.2e1 * t1076 + pkin(1)) * t1102 - t1148) * t1142 - qJ(3,2) * (t1070 * t1102 - t1101 * t1148)) * t1130 - (-t1100 * t1093 * t1161 - ((0.2e1 * t1075 + pkin(1)) * t1100 - t1149) * t1143 - qJ(3,3) * (t1069 * t1100 - t1099 * t1149)) * t1131) * m(3);];
taugX  = t1;
