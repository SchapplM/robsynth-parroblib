% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR1G1A0
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
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR1G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:00
% EndTime: 2020-03-09 20:34:02
% DurationCPUTime: 1.93s
% Computational Cost: add. (2118->281), mult. (5221->551), div. (2172->10), fcn. (4929->24), ass. (0->240)
t134 = sin(qJ(3,1));
t140 = cos(qJ(3,1));
t158 = rSges(3,1) * t140 - rSges(3,2) * t134;
t265 = m(3) * t158;
t132 = sin(qJ(3,2));
t138 = cos(qJ(3,2));
t159 = rSges(3,1) * t138 - rSges(3,2) * t132;
t264 = m(3) * t159;
t130 = sin(qJ(3,3));
t136 = cos(qJ(3,3));
t160 = rSges(3,1) * t136 - rSges(3,2) * t130;
t263 = m(3) * t160;
t104 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t262 = 0.2e1 * t104;
t151 = rSges(3,2) ^ 2;
t152 = rSges(3,1) ^ 2;
t95 = (-t151 + t152) * m(3) - Icges(3,1) + Icges(3,2);
t261 = t95 / 0.2e1;
t145 = m(2) * rSges(2,1);
t88 = t130 * rSges(3,1) + t136 * rSges(3,2);
t260 = m(3) * t88;
t89 = t132 * rSges(3,1) + t138 * rSges(3,2);
t259 = m(3) * t89;
t90 = t134 * rSges(3,1) + t140 * rSges(3,2);
t258 = m(3) * t90;
t257 = rSges(2,1) * g(3);
t256 = rSges(3,3) * m(3);
t255 = g(3) * rSges(3,3);
t102 = rSges(3,2) * t256 - Icges(3,6);
t254 = -t102 / 0.4e1;
t103 = rSges(3,1) * t256 - Icges(3,5);
t253 = t103 / 0.4e1;
t153 = 0.1e1 / pkin(2);
t252 = m(3) * t153;
t111 = m(1) + m(2) + m(3);
t251 = g(3) * t111;
t116 = 0.1e1 / t136;
t101 = m(2) * rSges(2,2) - t256;
t131 = sin(qJ(2,3));
t137 = cos(qJ(2,3));
t58 = (t145 + t263) * t137 - t131 * t101;
t249 = t116 * t58;
t124 = legFrame(3,3);
t105 = sin(t124);
t108 = cos(t124);
t146 = xDP(2);
t147 = xDP(1);
t221 = t116 * t153;
t64 = (t105 * t147 - t108 * t146) * t221;
t61 = t64 ^ 2;
t248 = t116 * t61;
t119 = 0.1e1 / t138;
t133 = sin(qJ(2,2));
t139 = cos(qJ(2,2));
t59 = (t145 + t264) * t139 - t133 * t101;
t247 = t119 * t59;
t125 = legFrame(2,3);
t106 = sin(t125);
t109 = cos(t125);
t219 = t119 * t153;
t65 = (t106 * t147 - t109 * t146) * t219;
t62 = t65 ^ 2;
t246 = t119 * t62;
t122 = 0.1e1 / t140;
t135 = sin(qJ(2,1));
t141 = cos(qJ(2,1));
t60 = (t145 + t265) * t141 - t135 * t101;
t245 = t122 * t60;
t126 = legFrame(1,3);
t107 = sin(t126);
t110 = cos(t126);
t217 = t122 * t153;
t66 = (t107 * t147 - t110 * t146) * t217;
t63 = t66 ^ 2;
t244 = t122 * t63;
t243 = t130 * t64;
t242 = t131 * t88;
t241 = t132 * t65;
t240 = t133 * t89;
t239 = t134 * t66;
t238 = t135 * t90;
t206 = t151 + t152;
t100 = m(3) * t206 + Icges(3,3);
t237 = t100 * t153;
t236 = t105 * t116;
t235 = t106 * t119;
t234 = t107 * t122;
t233 = t108 * t116;
t232 = t109 * t119;
t231 = t110 * t122;
t230 = t111 * t116;
t229 = t111 * t119;
t228 = t111 * t122;
t112 = 0.1e1 / t131;
t227 = t112 * t116;
t117 = 0.1e1 / t136 ^ 2;
t226 = t112 * t117;
t113 = 0.1e1 / t133;
t225 = t113 * t119;
t120 = 0.1e1 / t138 ^ 2;
t224 = t113 * t120;
t114 = 0.1e1 / t135;
t223 = t114 * t122;
t123 = 0.1e1 / t140 ^ 2;
t222 = t114 * t123;
t220 = t117 * t153;
t218 = t120 * t153;
t216 = t123 * t153;
t215 = t130 * t136;
t214 = t130 * t137;
t213 = t132 * t138;
t212 = t132 * t139;
t211 = t134 * t140;
t210 = t134 * t141;
t209 = t136 * t137;
t208 = t138 * t139;
t207 = t140 * t141;
t205 = 0.2e1 * m(3);
t189 = t116 * t242;
t169 = t189 * t252;
t67 = -t105 * t214 - t108 * t136;
t187 = t67 * t220;
t70 = t105 * t130 + t108 * t209;
t28 = -t105 * t169 + (t187 * t58 + t230 * t70) * t112;
t184 = t119 * t240;
t168 = t184 * t252;
t68 = -t106 * t212 - t109 * t138;
t182 = t68 * t218;
t72 = t106 * t132 + t109 * t208;
t29 = -t106 * t168 + (t182 * t59 + t229 * t72) * t113;
t179 = t122 * t238;
t167 = t179 * t252;
t69 = -t107 * t210 - t110 * t140;
t177 = t69 * t216;
t74 = t107 * t134 + t110 * t207;
t30 = -t107 * t167 + (t177 * t60 + t228 * t74) * t114;
t204 = t30 + t29 + t28;
t71 = -t136 * t105 + t108 * t214;
t186 = t71 * t220;
t76 = t105 * t209 - t108 * t130;
t31 = t108 * t169 + (t186 * t58 + t230 * t76) * t112;
t73 = -t138 * t106 + t109 * t212;
t181 = t73 * t218;
t77 = t106 * t208 - t109 * t132;
t32 = t109 * t168 + (t181 * t59 + t229 * t77) * t113;
t75 = -t140 * t107 + t110 * t210;
t176 = t75 * t216;
t78 = t107 * t207 - t110 * t134;
t33 = t110 * t167 + (t176 * t60 + t228 * t78) * t114;
t203 = t33 + t32 + t31;
t202 = t70 * t227;
t201 = t76 * t227;
t200 = t67 * t226;
t199 = t71 * t226;
t198 = t72 * t225;
t197 = t77 * t225;
t196 = t68 * t224;
t195 = t73 * t224;
t194 = t74 * t223;
t193 = t78 * t223;
t192 = t69 * t222;
t191 = t75 * t222;
t190 = t130 * t248;
t79 = -t102 * t136 - t103 * t130;
t188 = t79 * t221;
t185 = t132 * t246;
t80 = -t102 * t138 - t103 * t132;
t183 = t80 * t219;
t180 = t134 * t244;
t81 = -t102 * t140 - t103 * t134;
t178 = t81 * t217;
t175 = t112 * t220;
t174 = t113 * t218;
t173 = t114 * t216;
t172 = t79 * t175;
t171 = t80 * t174;
t170 = t81 * t173;
t166 = t130 * t61 * t189;
t165 = t132 * t62 * t184;
t164 = t134 * t63 * t179;
t85 = t108 * g(1) + t105 * g(2);
t163 = g(3) * t131 + t137 * t85;
t86 = t109 * g(1) + t106 * g(2);
t162 = g(3) * t133 + t139 * t86;
t87 = t110 * g(1) + t107 * g(2);
t161 = g(3) * t135 + t141 * t87;
t52 = (t146 * t71 + t147 * t67) * t175;
t13 = ((-t131 * t243 + t209 * t52) * t116 * t52 + (-t131 * t215 * t52 + t137 * t64) * t117 * t64) * t112;
t49 = t52 ^ 2;
t37 = (-t136 * t49 - t248) * t112 * pkin(2);
t157 = -t137 * (t205 * t64 * t88 + t52 * t101) * t52 - t58 * t13 + (-t49 * t145 - (t49 + t61) * t263) * t131 - t111 * t37;
t53 = (t146 * t73 + t147 * t68) * t174;
t14 = ((-t133 * t241 + t208 * t53) * t119 * t53 + (-t133 * t213 * t53 + t139 * t65) * t120 * t65) * t113;
t50 = t53 ^ 2;
t38 = (-t138 * t50 - t246) * t113 * pkin(2);
t156 = -t139 * (t205 * t65 * t89 + t53 * t101) * t53 - t59 * t14 + (-t50 * t145 - (t50 + t62) * t264) * t133 - t111 * t38;
t54 = (t146 * t75 + t147 * t69) * t173;
t15 = ((-t135 * t239 + t207 * t54) * t122 * t54 + (-t135 * t211 * t54 + t141 * t66) * t123 * t66) * t114;
t51 = t54 ^ 2;
t39 = (-t140 * t51 - t244) * t114 * pkin(2);
t155 = -t141 * (t205 * t66 * t90 + t54 * t101) * t54 - t60 * t15 + (-t51 * t145 - (t51 + t63) * t265) * t135 - t111 * t39;
t154 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (0.2e1 * rSges(3,3) ^ 2 + t206) * m(3) / 0.2e1;
t150 = 0.2e1 * qJ(3,1);
t149 = 0.2e1 * qJ(3,2);
t148 = 0.2e1 * qJ(3,3);
t144 = rSges(2,2) * g(3);
t129 = xDDP(1);
t128 = xDDP(2);
t127 = xDDP(3);
t121 = t140 ^ 2;
t118 = t138 ^ 2;
t115 = t136 ^ 2;
t84 = -t107 * g(1) + t110 * g(2);
t83 = -t106 * g(1) + t109 * g(2);
t82 = -t105 * g(1) + t108 * g(2);
t57 = cos(t150) * t261 - t104 * sin(t150) + t154;
t56 = cos(t149) * t261 - t104 * sin(t149) + t154;
t55 = cos(t148) * t261 - t104 * sin(t148) + t154;
t45 = t75 * t170 + (-t110 * t237 - t258 * t78) * t122;
t44 = t73 * t171 + (-t109 * t237 - t259 * t77) * t119;
t43 = t71 * t172 + (-t108 * t237 - t260 * t76) * t116;
t42 = t69 * t170 + (t107 * t237 - t258 * t74) * t122;
t41 = t68 * t171 + (t106 * t237 - t259 * t72) * t119;
t40 = t67 * t172 + (t105 * t237 - t260 * t70) * t116;
t24 = -t110 * t178 + (t176 * t57 + t245 * t78) * t114;
t23 = -t109 * t183 + (t181 * t56 + t247 * t77) * t113;
t22 = -t108 * t188 + (t186 * t55 + t249 * t76) * t112;
t21 = t107 * t178 + (t177 * t57 + t245 * t74) * t114;
t20 = t106 * t183 + (t182 * t56 + t247 * t72) * t113;
t19 = t105 * t188 + (t187 * t55 + t249 * t70) * t112;
t9 = -t81 * t15 + t100 * t180 + t51 * (t121 * t262 + t211 * t95 + Icges(3,4)) + (t39 * t238 + (-t134 * t84 + t140 * t161) * rSges(3,2) + (-t51 * rSges(3,2) + t134 * t161 + t84 * t140) * rSges(3,1)) * m(3);
t8 = -t80 * t14 + t100 * t185 + t50 * (t118 * t262 + t213 * t95 + Icges(3,4)) + (t38 * t240 + (-t132 * t83 + t138 * t162) * rSges(3,2) + (-t50 * rSges(3,2) + t132 * t162 + t83 * t138) * rSges(3,1)) * m(3);
t7 = -t79 * t13 + t100 * t190 + t49 * (t115 * t262 + t215 * t95 + Icges(3,4)) + (t37 * t242 + (-t130 * t82 + t136 * t163) * rSges(3,2) + (-t49 * rSges(3,2) + t130 * t163 + t82 * t136) * rSges(3,1)) * m(3);
t6 = -m(3) * t164 + t155 - t251;
t5 = -m(3) * t165 + t156 - t251;
t4 = -m(3) * t166 + t157 - t251;
t3 = -t60 * t39 - t57 * t15 + t81 * t180 - 0.4e1 * ((t134 * t261 * t54 + t253 * t66) * t140 + t239 * t254 + (t121 - 0.1e1 / 0.2e1) * t54 * t104) * t66 + (m(2) * (rSges(2,2) * t87 - t257) + (-t87 * rSges(3,3) - g(3) * t158) * m(3)) * t141 + (m(2) * (rSges(2,1) * t87 + t144) + (t158 * t87 - t255) * m(3)) * t135;
t2 = -t59 * t38 - t56 * t14 + t80 * t185 - 0.4e1 * ((t132 * t261 * t53 + t253 * t65) * t138 + t241 * t254 + (t118 - 0.1e1 / 0.2e1) * t53 * t104) * t65 + (m(2) * (rSges(2,2) * t86 - t257) + (-t86 * rSges(3,3) - g(3) * t159) * m(3)) * t139 + (m(2) * (rSges(2,1) * t86 + t144) + (t159 * t86 - t255) * m(3)) * t133;
t1 = -t58 * t37 - t55 * t13 + t79 * t190 - 0.4e1 * ((t130 * t261 * t52 + t253 * t64) * t136 + t243 * t254 + (t115 - 0.1e1 / 0.2e1) * t52 * t104) * t64 + (m(2) * (rSges(2,2) * t85 - t257) + (-t85 * rSges(3,3) - g(3) * t160) * m(3)) * t137 + (m(2) * (rSges(2,1) * t85 + t144) + (t160 * t85 - t255) * m(3)) * t131;
t10 = [t4 * t202 + t5 * t198 + t6 * t194 - m(4) * g(1) + (t30 * t193 + t29 * t197 + t28 * t201) * t128 + t204 * t127 + (t30 * t194 + t29 * t198 + t28 * t202 + m(4)) * t129 + (t1 * t200 + t2 * t196 + t3 * t192 + t7 * t236 + t8 * t235 + t9 * t234 + (t19 * t200 + t192 * t21 + t196 * t20 + t234 * t42 + t235 * t41 + t236 * t40) * t129 + (t19 * t199 + t191 * t21 + t195 * t20 - t231 * t42 - t232 * t41 - t233 * t40) * t128) * t153; t4 * t201 + t5 * t197 + t6 * t193 - m(4) * g(2) + (t33 * t194 + t32 * t198 + t31 * t202) * t129 + t203 * t127 + (t33 * t193 + t32 * t197 + t31 * t201 + m(4)) * t128 + (t1 * t199 + t2 * t195 + t3 * t191 - t7 * t233 - t8 * t232 - t9 * t231 + (t192 * t24 + t196 * t23 + t200 * t22 + t234 * t45 + t235 * t44 + t236 * t43) * t129 + (t191 * t24 + t195 * t23 + t199 * t22 - t231 * t45 - t232 * t44 - t233 * t43) * t128) * t153; t204 * t129 + t203 * t128 + (-t164 - t165 - t166) * m(3) + t155 + t156 + t157 + (t127 - g(3)) * (0.3e1 * t111 + m(4));];
tauX  = t10;
