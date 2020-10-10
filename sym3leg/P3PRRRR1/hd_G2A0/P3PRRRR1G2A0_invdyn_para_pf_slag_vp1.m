% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR1G2A0
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
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR1G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2A0_invdyn_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR1G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:21
% EndTime: 2020-03-09 21:16:23
% DurationCPUTime: 2.10s
% Computational Cost: add. (2502->298), mult. (6180->595), div. (3120->10), fcn. (5337->24), ass. (0->218)
t122 = sin(qJ(2,3));
t103 = 0.1e1 / t122;
t127 = cos(qJ(3,3));
t107 = 0.1e1 / t127;
t138 = xDP(3);
t146 = 0.1e1 / pkin(2);
t139 = xDP(2);
t140 = xDP(1);
t115 = legFrame(3,2);
t96 = sin(t115);
t99 = cos(t115);
t163 = t139 * t96 - t140 * t99;
t108 = 0.1e1 / t127 ^ 2;
t121 = sin(qJ(3,3));
t128 = cos(qJ(2,3));
t186 = t108 * t121 * t128;
t46 = (-t107 * t138 - t163 * t186) * t146 * t103;
t40 = t46 ^ 2;
t124 = sin(qJ(2,2));
t104 = 0.1e1 / t124;
t129 = cos(qJ(3,2));
t110 = 0.1e1 / t129;
t116 = legFrame(2,2);
t100 = cos(t116);
t97 = sin(t116);
t162 = t100 * t140 - t139 * t97;
t111 = 0.1e1 / t129 ^ 2;
t123 = sin(qJ(3,2));
t130 = cos(qJ(2,2));
t185 = t111 * t123 * t130;
t47 = (-t110 * t138 + t162 * t185) * t146 * t104;
t41 = t47 ^ 2;
t126 = sin(qJ(2,1));
t105 = 0.1e1 / t126;
t131 = cos(qJ(3,1));
t113 = 0.1e1 / t131;
t117 = legFrame(1,2);
t101 = cos(t117);
t98 = sin(t117);
t161 = t101 * t140 - t139 * t98;
t114 = 0.1e1 / t131 ^ 2;
t125 = sin(qJ(3,1));
t132 = cos(qJ(2,1));
t184 = t114 * t125 * t132;
t48 = (-t113 * t138 + t161 * t184) * t146 * t105;
t42 = t48 ^ 2;
t95 = (m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4));
t242 = 2 * t95;
t144 = rSges(3,2) ^ 2;
t145 = rSges(3,1) ^ 2;
t86 = (-t144 + t145) * m(3) - Icges(3,1) + Icges(3,2);
t241 = t86 / 0.2e1;
t235 = rSges(3,3) * m(3);
t93 = rSges(3,2) * t235 - Icges(3,6);
t240 = -t93 / 0.4e1;
t94 = rSges(3,1) * t235 - Icges(3,5);
t239 = t94 / 0.4e1;
t137 = m(2) * rSges(2,1);
t79 = t121 * rSges(3,1) + t127 * rSges(3,2);
t238 = m(3) * t79;
t80 = t123 * rSges(3,1) + t129 * rSges(3,2);
t237 = m(3) * t80;
t81 = t125 * rSges(3,1) + t131 * rSges(3,2);
t236 = m(3) * t81;
t234 = g(3) * rSges(3,3);
t233 = rSges(3,1) * t127;
t232 = rSges(3,1) * t129;
t231 = rSges(3,1) * t131;
t169 = -rSges(3,2) * t121 + t233;
t92 = m(2) * rSges(2,2) - t235;
t55 = (m(3) * t169 + t137) * t128 - t122 * t92;
t230 = t107 * t55;
t209 = t107 * t146;
t61 = t163 * t209;
t58 = t61 ^ 2;
t229 = t107 * t58;
t208 = t110 * t146;
t62 = t162 * t208;
t59 = t62 ^ 2;
t228 = t110 * t59;
t202 = t124 * t129;
t71 = t100 * t202 + t97 * t123;
t227 = t110 * t71;
t226 = t110 * t80;
t207 = t113 * t146;
t63 = t161 * t207;
t60 = t63 ^ 2;
t225 = t113 * t60;
t201 = t126 * t131;
t72 = t101 * t201 + t98 * t125;
t224 = t113 * t72;
t223 = t113 * t81;
t222 = t122 * t79;
t221 = t128 * t61;
t220 = t130 * t62;
t219 = t132 * t63;
t200 = t144 + t145;
t91 = t200 * m(3) + Icges(3,3);
t218 = t146 * t91;
t217 = t146 * t97;
t216 = t146 * t98;
t215 = t100 * t146;
t214 = t101 * t146;
t102 = m(1) + m(2) + m(3);
t213 = t102 * t107;
t212 = t103 * t107;
t211 = t104 * t110;
t210 = t105 * t113;
t118 = xDDP(3);
t206 = t118 * t128;
t205 = t118 * t130;
t204 = t118 * t132;
t203 = t122 * t127;
t199 = m(3) * t226;
t198 = m(3) * t223;
t68 = -t100 * t123 + t97 * t202;
t197 = t68 * t211;
t69 = -t101 * t125 + t98 * t201;
t196 = t69 * t210;
t195 = t121 * t229;
t64 = -t94 * t121 - t93 * t127;
t194 = t64 * t209;
t193 = t123 * t228;
t65 = -t94 * t123 - t93 * t129;
t192 = t65 * t208;
t191 = t125 * t225;
t66 = -t94 * t125 - t93 * t131;
t190 = t66 * t207;
t189 = t86 * t121 * t127;
t188 = t86 * t123 * t129;
t187 = t86 * t125 * t131;
t183 = t124 * t199;
t182 = t126 * t198;
t181 = t103 * t186;
t180 = t104 * t185;
t179 = t105 * t184;
t178 = t146 * t186;
t177 = m(3) * t209 * t222;
t73 = t96 * g(1) + t99 * g(2);
t176 = g(3) * t128 + t122 * t73;
t74 = t97 * g(1) + t100 * g(2);
t175 = g(3) * t130 + t124 * t74;
t75 = t98 * g(1) + t101 * g(2);
t174 = g(3) * t132 + t126 * t75;
t173 = t96 * t178;
t172 = t99 * t178;
t171 = t185 * t217;
t170 = t184 * t216;
t168 = -rSges(3,2) * t123 + t232;
t167 = -rSges(3,2) * t125 + t231;
t119 = xDDP(2);
t120 = xDDP(1);
t67 = -t99 * t121 + t96 * t203;
t70 = t96 * t121 + t99 * t203;
t166 = t119 * t70 + t120 * t67;
t165 = t119 * t71 + t120 * t68;
t164 = t119 * t72 + t120 * t69;
t160 = t103 * t64 * t178;
t10 = ((-t61 * t122 * t121 + t127 * t46 * t128) * t107 * t46 + (-t121 * t46 * t203 + t221) * t108 * t61) * t103;
t106 = t127 ^ 2;
t135 = rSges(2,2) * g(3);
t136 = rSges(2,1) * g(3);
t28 = (-t127 * t40 - t229) * t103 * pkin(2);
t141 = 0.2e1 * qJ(3,3);
t147 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + ((2 * rSges(3,3) ^ 2) + t200) * m(3) / 0.2e1;
t43 = cos(t141) * t241 - t95 * sin(t141) + t147;
t1 = -t55 * t28 - t43 * t10 + t64 * t195 - 0.4e1 * ((t121 * t240 + t127 * t239) * t61 + (t189 / 0.2e1 + (t106 - 0.1e1 / 0.2e1) * t95) * t46) * t61 + (m(2) * (-rSges(2,1) * t73 + t135) + (-t169 * t73 - t234) * m(3)) * t128 + (m(2) * (rSges(2,2) * t73 + t136) + (-t73 * rSges(3,3) + g(3) * t169) * m(3)) * t122;
t76 = t99 * g(1) - t96 * g(2);
t159 = t1 * t181 - t107 * (-t64 * t10 + t91 * t195 + t40 * (t106 * t242 + Icges(3,4) + t189) + (t28 * t222 + (-t121 * t76 + t127 * t176) * rSges(3,2) + (-t40 * rSges(3,2) + t121 * t176 + t76 * t127) * rSges(3,1)) * m(3));
t109 = t129 ^ 2;
t11 = ((t62 * t124 * t123 + t129 * t47 * t130) * t110 * t47 - (-t123 * t47 * t202 - t220) * t111 * t62) * t104;
t29 = (-t129 * t41 - t228) * t104 * pkin(2);
t142 = 0.2e1 * qJ(3,2);
t44 = cos(t142) * t241 - t95 * sin(t142) + t147;
t56 = (m(3) * t168 + t137) * t130 - t124 * t92;
t2 = -t56 * t29 - t44 * t11 + t65 * t193 + 0.4e1 * (-(t123 * t240 + t129 * t239) * t62 + (t188 / 0.2e1 + (t109 - 0.1e1 / 0.2e1) * t95) * t47) * t62 + (m(2) * (-rSges(2,1) * t74 + t135) + (-t168 * t74 - t234) * m(3)) * t130 + (m(2) * (rSges(2,2) * t74 + t136) + (-t74 * rSges(3,3) + g(3) * t168) * m(3)) * t124;
t77 = t100 * g(1) - t97 * g(2);
t158 = -t110 * (-t65 * t11 + t91 * t193 + t41 * (t109 * t242 + Icges(3,4) + t188) + (t124 * t80 * t29 + (-t123 * t77 + t129 * t175) * rSges(3,2) + (-t41 * rSges(3,2) + t123 * t175 + t77 * t129) * rSges(3,1)) * m(3)) + t180 * t2;
t112 = t131 ^ 2;
t12 = ((t63 * t126 * t125 + t131 * t48 * t132) * t113 * t48 - (-t125 * t48 * t201 - t219) * t114 * t63) * t105;
t30 = (-t131 * t42 - t225) * t105 * pkin(2);
t143 = 0.2e1 * qJ(3,1);
t45 = cos(t143) * t241 - t95 * sin(t143) + t147;
t57 = (m(3) * t167 + t137) * t132 - t126 * t92;
t3 = -t57 * t30 - t45 * t12 + t66 * t191 + 0.4e1 * (-(t125 * t240 + t131 * t239) * t63 + (t187 / 0.2e1 + (t112 - 0.1e1 / 0.2e1) * t95) * t48) * t63 + (m(2) * (-rSges(2,1) * t75 + t135) + (-t167 * t75 - t234) * m(3)) * t132 + (m(2) * (rSges(2,2) * t75 + t136) + (-t75 * rSges(3,3) + g(3) * t167) * m(3)) * t126;
t78 = t101 * g(1) - t98 * g(2);
t157 = -t113 * (-t66 * t12 + t91 * t191 + t42 * (t112 * t242 + Icges(3,4) + t187) + (t126 * t81 * t30 + (-t125 * t78 + t131 * t174) * rSges(3,2) + (-t42 * rSges(3,2) + t125 * t174 + t78 * t131) * rSges(3,1)) * m(3)) + t179 * t3;
t13 = -t99 * t194 + (t172 * t43 + t67 * t230) * t103;
t156 = -t107 * (t99 * t160 + (-t99 * t218 - t67 * t238) * t107) + t13 * t181;
t16 = t96 * t194 + (-t173 * t43 + t70 * t230) * t103;
t155 = -t107 * (-t96 * t160 + (t96 * t218 - t70 * t238) * t107) + t16 * t181;
t25 = (t128 * t55 - t43 * t209) * t103;
t154 = -t107 * (-t103 * t194 - t128 * t238) + t181 * t25;
t14 = t56 * t197 + (-t110 * t65 + t180 * t44) * t215;
t153 = -t110 * (-t68 * t199 + (-t110 * t91 + t180 * t65) * t215) + t14 * t180;
t17 = t97 * t192 + (-t171 * t44 + t56 * t227) * t104;
t152 = -t110 * (-t104 * t65 * t171 + (t91 * t217 - t71 * t237) * t110) + t17 * t180;
t26 = (t130 * t56 - t44 * t208) * t104;
t151 = -t110 * (-t104 * t192 - t130 * t237) + t180 * t26;
t15 = t57 * t196 + (-t113 * t66 + t179 * t45) * t214;
t150 = -t113 * (-t69 * t198 + (-t113 * t91 + t179 * t66) * t214) + t15 * t179;
t18 = t98 * t190 + (-t170 * t45 + t57 * t224) * t105;
t149 = -t113 * (-t105 * t66 * t170 + (t91 * t216 - t72 * t236) * t113) + t179 * t18;
t27 = (t132 * t57 - t45 * t207) * t105;
t148 = -t113 * (-t105 * t190 - t132 * t236) + t179 * t27;
t39 = t42 + t60;
t38 = t41 + t59;
t37 = t40 + t58;
t24 = -t182 * t216 + (t102 * t224 - t170 * t57) * t105;
t23 = -t183 * t217 + (t102 * t227 - t171 * t56) * t104;
t22 = -t96 * t177 + (-t173 * t55 + t70 * t213) * t103;
t21 = t102 * t196 + (t179 * t57 + t182) * t214;
t20 = t102 * t197 + (t180 * t56 + t183) * t215;
t19 = t99 * t177 + (t172 * t55 + t67 * t213) * t103;
t6 = -t57 * t12 + (-t30 - t75) * t102 + (0.2e1 * t48 * t81 * t219 + (-t39 * t231 + (rSges(3,2) * t39 - t60 * t223) * t125) * t126) * m(3) + (-t126 * t137 - t132 * t92) * t42;
t5 = -t56 * t11 + (-t29 - t74) * t102 + (0.2e1 * t47 * t80 * t220 + (-t38 * t232 + (rSges(3,2) * t38 - t59 * t226) * t123) * t124) * m(3) + (-t124 * t137 - t130 * t92) * t41;
t4 = -t55 * t10 + (-t28 - t73) * t102 + (-0.2e1 * t46 * t79 * t221 + (-t37 * t233 + (rSges(3,2) * t37 - t79 * t229) * t121) * t122) * m(3) + (-t122 * t137 - t128 * t92) * t40;
t7 = [(-g(1) + t120) * m(4) + (t21 * t204 + (t164 * t21 + t6 * t69) * t113) * t105 + (t20 * t205 + (t165 * t20 + t5 * t68) * t110) * t104 + (t19 * t206 + (t166 * t19 + t4 * t67) * t107) * t103 + ((-t150 * t98 - t153 * t97 - t156 * t96) * t119 + (-t13 * t212 - t14 * t211 - t15 * t210) * t118 + (t120 * t156 + t159) * t99 + (t120 * t150 + t157) * t101 + (t120 * t153 + t158) * t100) * t146; (-g(2) + t119) * m(4) + (t24 * t204 + (t164 * t24 + t6 * t72) * t113) * t105 + (t23 * t205 + (t165 * t23 + t5 * t71) * t110) * t104 + (t22 * t206 + (t166 * t22 + t4 * t70) * t107) * t103 + ((t100 * t152 + t101 * t149 + t155 * t99) * t120 + (-t16 * t212 - t17 * t211 - t18 * t210) * t118 + (-t119 * t149 - t157) * t98 + (-t119 * t152 - t158) * t97 + (-t119 * t155 - t159) * t96) * t146; (-g(3) + t118) * m(4) + (t132 * t6 + (t113 * t164 + t204) * (t102 * t132 - t57 * t207) * t105) * t105 + (t130 * t5 + (t110 * t165 + t205) * (t102 * t130 - t56 * t208) * t104) * t104 + (t128 * t4 + (t107 * t166 + t206) * (t102 * t128 - t55 * t209) * t103) * t103 + (-t1 * t212 - t2 * t211 - t3 * t210 + (t100 * t151 + t101 * t148 + t154 * t99) * t120 + (-t148 * t98 - t151 * t97 - t154 * t96) * t119 + (-t27 * t210 - t26 * t211 - t25 * t212) * t118) * t146;];
tauX  = t7;
