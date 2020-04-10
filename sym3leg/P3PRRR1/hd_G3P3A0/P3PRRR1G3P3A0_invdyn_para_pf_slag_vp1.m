% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRR1G3P3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:42
% EndTime: 2020-03-09 21:06:45
% DurationCPUTime: 2.30s
% Computational Cost: add. (17271->269), mult. (9210->498), div. (3366->5), fcn. (9498->18), ass. (0->190)
t142 = 0.1e1 / pkin(2);
t124 = pkin(7) + qJ(2,1);
t113 = qJ(3,1) + t124;
t100 = sin(t113);
t103 = cos(t113);
t107 = sin(t124);
t110 = cos(t124);
t220 = 0.1e1 / (t100 * t110 - t103 * t107);
t225 = t142 * t220;
t123 = pkin(7) + qJ(2,2);
t112 = qJ(3,2) + t123;
t102 = cos(t112);
t106 = sin(t123);
t109 = cos(t123);
t99 = sin(t112);
t221 = 0.1e1 / (-t102 * t106 + t109 * t99);
t224 = t142 * t221;
t122 = pkin(7) + qJ(2,3);
t111 = qJ(3,3) + t122;
t101 = cos(t111);
t105 = sin(t122);
t108 = cos(t122);
t98 = sin(t111);
t222 = 0.1e1 / (-t101 * t105 + t108 * t98);
t223 = t142 * t222;
t219 = rSges(3,1) * g(3);
t218 = rSges(2,2) * g(3);
t217 = rSges(3,2) * g(3);
t134 = xDP(3);
t140 = 0.1e1 / pkin(3);
t136 = xDP(1);
t127 = legFrame(3,2);
t117 = cos(t127);
t70 = pkin(2) * t108 + pkin(3) * t101;
t180 = t117 * t222 * t70;
t164 = t136 * t180;
t160 = t140 * t164;
t114 = sin(t127);
t135 = xDP(2);
t193 = t114 * t135;
t67 = pkin(2) * t105 + pkin(3) * t98;
t200 = t67 * t140;
t210 = t140 * t70;
t183 = t136 * t142;
t197 = t101 * t117;
t49 = t222 * t183 * t197;
t16 = t49 + (-t160 + ((-t98 + t200) * t134 + (-t101 + t210) * t193) * t222) * t142;
t182 = t140 * t142;
t34 = (-t164 + (t134 * t67 + t70 * t193) * t222) * t182;
t216 = t16 * t34;
t128 = legFrame(2,2);
t118 = cos(t128);
t71 = pkin(2) * t109 + pkin(3) * t102;
t179 = t118 * t221 * t71;
t163 = t136 * t179;
t159 = t140 * t163;
t115 = sin(t128);
t191 = t115 * t135;
t68 = pkin(2) * t106 + pkin(3) * t99;
t199 = t68 * t140;
t209 = t140 * t71;
t196 = t102 * t118;
t50 = t221 * t183 * t196;
t17 = t50 + (-t159 + ((-t99 + t199) * t134 + (-t102 + t209) * t191) * t221) * t142;
t35 = (-t163 + (t134 * t68 + t71 * t191) * t221) * t182;
t215 = t17 * t35;
t129 = legFrame(1,2);
t119 = cos(t129);
t72 = pkin(2) * t110 + pkin(3) * t103;
t178 = t119 * t220 * t72;
t162 = t136 * t178;
t158 = t140 * t162;
t116 = sin(t129);
t189 = t116 * t135;
t69 = pkin(2) * t107 + pkin(3) * t100;
t198 = t69 * t140;
t208 = t140 * t72;
t195 = t103 * t119;
t51 = t220 * t183 * t195;
t18 = t51 + (-t158 + ((-t100 + t198) * t134 + (-t103 + t208) * t189) * t220) * t142;
t36 = (-t162 + (t69 * t134 + t72 * t189) * t220) * t182;
t214 = t18 * t36;
t213 = t101 * t222;
t212 = t102 * t221;
t211 = t103 * t220;
t120 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t88 = m(3) * t120 + Icges(3,3);
t207 = t140 * t88;
t131 = xDDP(2);
t194 = t114 * t131;
t192 = t115 * t131;
t190 = t116 * t131;
t188 = t117 * t142;
t187 = t118 * t142;
t186 = t119 * t142;
t185 = t131 * t140;
t132 = xDDP(1);
t184 = t132 * t140;
t125 = m(1) + m(2) + m(3);
t181 = (t114 * t117 + t115 * t118 + t116 * t119) * t125;
t73 = -rSges(3,1) * t108 - t105 * rSges(3,2);
t74 = rSges(3,1) * t105 - rSges(3,2) * t108;
t148 = (-t101 * t73 + t74 * t98) * pkin(2);
t46 = Icges(3,3) + (t120 + t148) * m(3);
t177 = t46 * t210;
t75 = -rSges(3,1) * t109 - t106 * rSges(3,2);
t76 = rSges(3,1) * t106 - rSges(3,2) * t109;
t147 = (-t102 * t75 + t76 * t99) * pkin(2);
t47 = Icges(3,3) + (t120 + t147) * m(3);
t176 = t47 * t209;
t77 = -rSges(3,1) * t110 - t107 * rSges(3,2);
t78 = rSges(3,1) * t107 - rSges(3,2) * t110;
t146 = (t100 * t78 - t103 * t77) * pkin(2);
t48 = Icges(3,3) + (t120 + t146) * m(3);
t175 = t48 * t208;
t174 = t70 * t207;
t173 = t71 * t207;
t172 = t72 * t207;
t171 = t114 * t223;
t170 = t115 * t224;
t169 = t116 * t225;
t168 = 0.2e1 * pkin(2) * pkin(3);
t167 = t132 * t197;
t166 = t132 * t196;
t165 = t132 * t195;
t161 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3);
t85 = g(1) * t117 - g(2) * t114;
t157 = t101 * (rSges(3,2) * t85 + t219) + (rSges(3,1) * t85 - t217) * t98;
t156 = t101 * t74 + t73 * t98;
t86 = g(1) * t118 - g(2) * t115;
t155 = t102 * (rSges(3,2) * t86 + t219) + (rSges(3,1) * t86 - t217) * t99;
t154 = t102 * t76 + t75 * t99;
t87 = g(1) * t119 - g(2) * t116;
t153 = (rSges(3,1) * t87 - t217) * t100 + t103 * (rSges(3,2) * t87 + t219);
t152 = t100 * t77 + t103 * t78;
t151 = t101 * t108 + t105 * t98;
t150 = t102 * t109 + t106 * t99;
t149 = t100 * t107 + t103 * t110;
t145 = t151 * pkin(2);
t144 = t150 * pkin(2);
t143 = t149 * pkin(2);
t141 = pkin(2) ^ 2;
t139 = pkin(3) ^ 2;
t133 = rSges(2,1) * g(3);
t130 = xDDP(3);
t126 = g(3) * pkin(2) * m(3);
t104 = t141 + t120;
t84 = g(1) * t116 + g(2) * t119;
t83 = g(1) * t115 + g(2) * t118;
t82 = g(1) * t114 + g(2) * t117;
t45 = (t104 + 0.2e1 * t146) * m(3) + t161;
t44 = (t104 + 0.2e1 * t147) * m(3) + t161;
t43 = (t104 + 0.2e1 * t148) * m(3) + t161;
t42 = t51 + (-t100 * t134 - t103 * t189) * t225;
t41 = t50 + (-t102 * t191 - t99 * t134) * t224;
t40 = t49 + (-t101 * t193 - t98 * t134) * t223;
t39 = (-t100 * t48 + t88 * t198) * t225;
t38 = (t88 * t199 - t47 * t99) * t224;
t37 = (t88 * t200 - t46 * t98) * t223;
t33 = (-t172 * t220 + t48 * t211) * t186;
t32 = (-t173 * t221 + t47 * t212) * t187;
t31 = (-t174 * t222 + t46 * t213) * t188;
t30 = (-t103 * t48 + t172) * t169;
t29 = (-t102 * t47 + t173) * t170;
t28 = (-t101 * t46 + t174) * t171;
t27 = (-t100 * t45 + t48 * t198) * t225;
t26 = (t47 * t199 - t44 * t99) * t224;
t25 = (t46 * t200 - t43 * t98) * t223;
t24 = (-t175 * t220 + t45 * t211) * t186;
t23 = (-t176 * t221 + t44 * t212) * t187;
t22 = (-t177 * t222 + t43 * t213) * t188;
t21 = (-t103 * t45 + t175) * t169;
t20 = (-t102 * t44 + t176) * t170;
t19 = (-t101 * t43 + t177) * t171;
t15 = t51 + (-t158 / 0.2e1 + ((-t100 + t198 / 0.2e1) * t134 + (-t103 + t208 / 0.2e1) * t189) * t220) * t142;
t14 = t50 + (-t159 / 0.2e1 + ((-t99 + t199 / 0.2e1) * t134 + (-t102 + t209 / 0.2e1) * t191) * t221) * t142;
t13 = t49 + (-t160 / 0.2e1 + ((-t98 + t200 / 0.2e1) * t134 + (-t101 + t210 / 0.2e1) * t193) * t222) * t142;
t12 = (pkin(3) * t214 + (t18 * pkin(3) + t42 * t143) * t42) * t225;
t11 = (pkin(3) * t215 + (t17 * pkin(3) + t41 * t144) * t41) * t224;
t10 = (pkin(3) * t216 + (t16 * pkin(3) + t40 * t145) * t40) * t223;
t9 = ((-t149 * t15 * t168 - t139 * t18 - t141 * t42) * t140 * t42 - (pkin(3) + t143) * t214) * t225;
t8 = ((-t150 * t14 * t168 - t139 * t17 - t141 * t41) * t140 * t41 - (pkin(3) + t144) * t215) * t224;
t7 = ((-t151 * t13 * t168 - t139 * t16 - t141 * t40) * t140 * t40 - (pkin(3) + t145) * t216) * t223;
t6 = t48 * t12 + t88 * t9 + (-t42 ^ 2 * t152 * pkin(2) + t153) * m(3);
t5 = t47 * t11 + t88 * t8 + (-t41 ^ 2 * t154 * pkin(2) + t155) * m(3);
t4 = t46 * t10 + t88 * t7 + (-t40 ^ 2 * t156 * pkin(2) + t157) * m(3);
t3 = t45 * t12 + t48 * t9 + (t126 + m(2) * (rSges(2,2) * t87 + t133)) * t110 + m(2) * (rSges(2,1) * t87 - t218) * t107 + ((0.2e1 * t152 * t36 * t15 + t87 * t107) * pkin(2) + t153) * m(3);
t2 = t44 * t11 + t47 * t8 + (t126 + m(2) * (rSges(2,2) * t86 + t133)) * t109 + m(2) * (rSges(2,1) * t86 - t218) * t106 + ((0.2e1 * t154 * t35 * t14 + t86 * t106) * pkin(2) + t155) * m(3);
t1 = t43 * t10 + t46 * t7 + (t126 + m(2) * (rSges(2,2) * t85 + t133)) * t108 + m(2) * (rSges(2,1) * t85 - t218) * t105 + ((0.2e1 * t156 * t34 * t13 + t85 * t105) * pkin(2) + t157) * m(3);
t52 = [(-g(1) + t132) * m(4) + t181 * t131 + (-t114 * t82 - t115 * t83 - t116 * t84 + (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) * t132) * t125 + ((-t100 * t130 * t24 + (-t24 * t190 + (t132 * t24 + t3) * t119) * t103) * t220 + (-t4 * t180 - t5 * t179 - t6 * t178 + (t130 * t69 + t72 * t190) * t220 * t33 + (-t33 * t178 - t32 * t179 - t31 * t180) * t132) * t140 + (t115 * t32 * t71 * t185 + (t32 * t199 - t23 * t99) * t130 + (-t23 * t192 + (t132 * t23 + t2) * t118) * t102) * t221 + (t114 * t31 * t70 * t185 + (t31 * t200 - t22 * t98) * t130 + (-t22 * t194 + (t132 * t22 + t1) * t117) * t101) * t222) * t142; (-g(2) + t131) * m(4) + t181 * t132 + (-t117 * t82 - t118 * t83 - t119 * t84 + (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) * t131) * t125 + ((-t30 * t178 - t29 * t179 - t28 * t180) * t184 + (t21 * t165 + (-t100 * t21 + t30 * t198) * t130 + ((-t103 * t21 + t30 * t208) * t131 - t103 * t3 + t6 * t208) * t116) * t220 + (t20 * t166 + (t29 * t199 - t20 * t99) * t130 + ((-t102 * t20 + t29 * t209) * t131 - t102 * t2 + t5 * t209) * t115) * t221 + (t19 * t167 + (-t19 * t98 + t28 * t200) * t130 + ((-t101 * t19 + t28 * t210) * t131 - t101 * t1 + t4 * t210) * t114) * t222) * t142; (-g(3) + t130) * m(4) + ((-t39 * t178 - t38 * t179 - t37 * t180) * t184 + (t27 * t165 + (-t100 * t27 + t39 * t198) * t130 - t100 * t3 + t6 * t198 + (-t103 * t27 + t39 * t208) * t190) * t220 + (t26 * t166 + (t38 * t199 - t26 * t99) * t130 - t99 * t2 + t5 * t199 + (-t102 * t26 + t38 * t209) * t192) * t221 + (t25 * t167 + (t37 * t200 - t25 * t98) * t130 - t98 * t1 + t4 * t200 + (-t101 * t25 + t37 * t210) * t194) * t222) * t142;];
tauX  = t52;
