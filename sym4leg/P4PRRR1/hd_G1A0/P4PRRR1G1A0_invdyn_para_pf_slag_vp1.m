% Calculate vector of inverse dynamics forces for parallel robot
% P4PRRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% xDDP [4x1]
%   Generalized platform accelerations
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tauX [4x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1P1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:10:23
% EndTime: 2020-03-02 20:10:30
% DurationCPUTime: 7.53s
% Computational Cost: add. (23591->342), mult. (15967->640), div. (3184->6), fcn. (16030->26), ass. (0->262)
t238 = 0.1e1 / pkin(2);
t208 = pkin(7) + qJ(2,1);
t192 = qJ(3,1) + t208;
t176 = sin(t192);
t179 = cos(t192);
t186 = sin(t208);
t189 = cos(t208);
t320 = 0.1e1 / (t176 * t189 - t179 * t186);
t327 = t320 * t238;
t207 = pkin(7) + qJ(2,2);
t191 = qJ(3,2) + t207;
t175 = sin(t191);
t178 = cos(t191);
t185 = sin(t207);
t188 = cos(t207);
t321 = 0.1e1 / (t175 * t188 - t178 * t185);
t326 = t321 * t238;
t206 = pkin(7) + qJ(2,3);
t190 = qJ(3,3) + t206;
t174 = sin(t190);
t177 = cos(t190);
t184 = sin(t206);
t187 = cos(t206);
t322 = 0.1e1 / (t174 * t187 - t177 * t184);
t325 = t322 * t238;
t205 = pkin(7) + qJ(2,4);
t183 = qJ(3,4) + t205;
t172 = sin(t183);
t173 = cos(t183);
t181 = sin(t205);
t182 = cos(t205);
t323 = 0.1e1 / (t172 * t182 - t173 * t181);
t324 = t323 * t238;
t236 = 0.1e1 / pkin(3);
t285 = t236 * t238;
t274 = t323 * t285;
t222 = xP(4);
t203 = sin(t222);
t204 = cos(t222);
t225 = koppelP(4,2);
t229 = koppelP(4,1);
t143 = t203 * t229 + t204 * t225;
t219 = xDP(4);
t221 = xDP(1);
t129 = -t143 * t219 + t221;
t211 = legFrame(4,3);
t193 = sin(t211);
t197 = cos(t211);
t118 = -t193 * t172 + t197 * t173;
t82 = -pkin(2) * (t181 * t193 - t182 * t197) + t118 * pkin(3);
t313 = t82 * t129;
t147 = -t203 * t225 + t204 * t229;
t220 = xDP(2);
t125 = t147 * t219 + t220;
t117 = t197 * t172 + t193 * t173;
t81 = pkin(2) * (t181 * t197 + t182 * t193) + t117 * pkin(3);
t314 = t81 * t125;
t242 = (t313 + t314) * t274;
t57 = (t117 * t125 + t118 * t129) * t324;
t26 = -t242 + t57;
t319 = t26 * t242;
t272 = t322 * t285;
t226 = koppelP(3,2);
t230 = koppelP(3,1);
t144 = t203 * t230 + t204 * t226;
t130 = -t144 * t219 + t221;
t212 = legFrame(3,3);
t194 = sin(t212);
t198 = cos(t212);
t120 = -t194 * t174 + t198 * t177;
t86 = -pkin(2) * (t184 * t194 - t187 * t198) + t120 * pkin(3);
t309 = t86 * t130;
t148 = -t203 * t226 + t204 * t230;
t126 = t148 * t219 + t220;
t119 = t198 * t174 + t194 * t177;
t83 = pkin(2) * (t184 * t198 + t187 * t194) + t119 * pkin(3);
t312 = t83 * t126;
t241 = (t309 + t312) * t272;
t58 = (t119 * t126 + t120 * t130) * t325;
t30 = -t241 + t58;
t318 = t30 * t241;
t270 = t321 * t285;
t227 = koppelP(2,2);
t231 = koppelP(2,1);
t145 = t203 * t231 + t204 * t227;
t131 = -t145 * t219 + t221;
t213 = legFrame(2,3);
t195 = sin(t213);
t199 = cos(t213);
t122 = -t195 * t175 + t199 * t178;
t87 = -pkin(2) * (t185 * t195 - t188 * t199) + t122 * pkin(3);
t308 = t87 * t131;
t149 = -t203 * t227 + t204 * t231;
t127 = t149 * t219 + t220;
t121 = t199 * t175 + t195 * t178;
t84 = pkin(2) * (t185 * t199 + t188 * t195) + t121 * pkin(3);
t311 = t84 * t127;
t240 = (t308 + t311) * t270;
t59 = (t121 * t127 + t122 * t131) * t326;
t31 = -t240 + t59;
t317 = t31 * t240;
t268 = t320 * t285;
t228 = koppelP(1,2);
t232 = koppelP(1,1);
t146 = t203 * t232 + t204 * t228;
t132 = -t146 * t219 + t221;
t214 = legFrame(1,3);
t196 = sin(t214);
t200 = cos(t214);
t124 = -t196 * t176 + t200 * t179;
t88 = -pkin(2) * (t186 * t196 - t189 * t200) + t124 * pkin(3);
t307 = t88 * t132;
t150 = -t203 * t228 + t204 * t232;
t128 = t150 * t219 + t220;
t123 = t200 * t176 + t196 * t179;
t85 = pkin(2) * (t186 * t200 + t189 * t196) + t123 * pkin(3);
t310 = t85 * t128;
t239 = (t307 + t310) * t268;
t60 = (t123 * t128 + t124 * t132) * t327;
t32 = -t239 + t60;
t316 = t32 * t239;
t306 = t323 * t117;
t305 = t323 * t118;
t303 = t323 * t236;
t301 = t322 * t119;
t300 = t322 * t120;
t298 = t321 * t121;
t297 = t321 * t122;
t295 = t320 * t123;
t294 = t320 * t124;
t292 = t322 * t236;
t290 = t321 * t236;
t288 = t320 * t236;
t201 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t159 = m(3) * t201 + Icges(3,3);
t286 = t159 * t236;
t284 = t81 * t303;
t283 = t82 * t303;
t282 = t83 * t292;
t281 = t86 * t292;
t280 = t84 * t290;
t279 = t87 * t290;
t278 = t85 * t288;
t277 = t88 * t288;
t276 = 0.2e1 * pkin(2) * pkin(3);
t275 = t323 * t286;
t273 = t322 * t286;
t271 = t321 * t286;
t269 = t320 * t286;
t267 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3);
t151 = -g(1) * t193 + g(2) * t197;
t155 = g(1) * t197 + g(2) * t193;
t266 = (rSges(3,1) * t155 + rSges(3,2) * t151) * t172 - t173 * (rSges(3,1) * t151 - rSges(3,2) * t155);
t152 = -g(1) * t194 + g(2) * t198;
t156 = g(1) * t198 + g(2) * t194;
t265 = (rSges(3,1) * t156 + rSges(3,2) * t152) * t174 - t177 * (rSges(3,1) * t152 - rSges(3,2) * t156);
t153 = -g(1) * t195 + g(2) * t199;
t157 = g(1) * t199 + g(2) * t195;
t264 = (rSges(3,1) * t157 + rSges(3,2) * t153) * t175 - t178 * (rSges(3,1) * t153 - rSges(3,2) * t157);
t154 = -g(1) * t196 + g(2) * t200;
t158 = g(1) * t200 + g(2) * t196;
t263 = -(rSges(3,1) * t154 - rSges(3,2) * t158) * t179 + (rSges(3,1) * t158 + rSges(3,2) * t154) * t176;
t133 = -rSges(3,1) * t182 - t181 * rSges(3,2);
t134 = rSges(3,1) * t181 - rSges(3,2) * t182;
t262 = t133 * t172 + t134 * t173;
t137 = -rSges(3,1) * t187 - t184 * rSges(3,2);
t138 = rSges(3,1) * t184 - rSges(3,2) * t187;
t261 = t137 * t174 + t138 * t177;
t139 = -rSges(3,1) * t188 - t185 * rSges(3,2);
t140 = rSges(3,1) * t185 - rSges(3,2) * t188;
t260 = t139 * t175 + t140 * t178;
t141 = -rSges(3,1) * t189 - t186 * rSges(3,2);
t142 = rSges(3,1) * t186 - rSges(3,2) * t189;
t259 = t141 * t176 + t142 * t179;
t258 = t172 * t181 + t173 * t182;
t257 = t174 * t184 + t177 * t187;
t256 = t175 * t185 + t178 * t188;
t255 = t176 * t186 + t179 * t189;
t210 = t219 ^ 2;
t215 = xDDP(4);
t217 = xDDP(2);
t89 = -t143 * t210 + t147 * t215 + t217;
t218 = xDDP(1);
t93 = -t143 * t215 - t147 * t210 + t218;
t254 = t238 * (t81 * t89 + t82 * t93);
t90 = -t144 * t210 + t148 * t215 + t217;
t94 = -t144 * t215 - t148 * t210 + t218;
t253 = t238 * (t83 * t90 + t86 * t94);
t91 = -t145 * t210 + t149 * t215 + t217;
t95 = -t145 * t215 - t149 * t210 + t218;
t252 = t238 * (t84 * t91 + t87 * t95);
t92 = -t146 * t210 + t150 * t215 + t217;
t96 = -t146 * t215 - t150 * t210 + t218;
t251 = t238 * (t85 * t92 + t88 * t96);
t250 = (-t133 * t173 + t134 * t172) * pkin(2);
t249 = (-t137 * t177 + t138 * t174) * pkin(2);
t248 = (-t139 * t178 + t140 * t175) * pkin(2);
t247 = (-t141 * t179 + t142 * t176) * pkin(2);
t246 = t258 * pkin(2);
t245 = t257 * pkin(2);
t244 = t256 * pkin(2);
t243 = t255 * pkin(2);
t237 = pkin(2) ^ 2;
t235 = pkin(3) ^ 2;
t224 = rSges(4,1);
t223 = rSges(4,2);
t180 = t237 + t201;
t136 = -t223 * t203 + t224 * t204;
t135 = t224 * t203 + t223 * t204;
t80 = Icges(3,3) + (t201 + t247) * m(3);
t79 = Icges(3,3) + (t201 + t248) * m(3);
t78 = Icges(3,3) + (t201 + t249) * m(3);
t77 = Icges(3,3) + (t201 + t250) * m(3);
t76 = (t180 + 0.2e1 * t247) * m(3) + t267;
t75 = (t180 + 0.2e1 * t248) * m(3) + t267;
t74 = (t180 + 0.2e1 * t249) * m(3) + t267;
t73 = (t180 + 0.2e1 * t250) * m(3) + t267;
t64 = (t123 * t150 - t124 * t146) * t327;
t63 = (t121 * t149 - t122 * t145) * t326;
t62 = (t119 * t148 - t120 * t144) * t325;
t61 = (t117 * t147 - t118 * t143) * t324;
t56 = (-t146 * t88 + t150 * t85) * t268;
t55 = (-t145 * t87 + t149 * t84) * t270;
t54 = (-t144 * t86 + t148 * t83) * t272;
t53 = (-t143 * t82 + t147 * t81) * t274;
t40 = (-t80 * t277 + t76 * t294) * t238;
t39 = (-t80 * t278 + t76 * t295) * t238;
t38 = (-t79 * t279 + t75 * t297) * t238;
t37 = (-t79 * t280 + t75 * t298) * t238;
t36 = (-t78 * t281 + t74 * t300) * t238;
t35 = (-t78 * t282 + t74 * t301) * t238;
t34 = (-t77 * t283 + t73 * t305) * t238;
t33 = (-t77 * t284 + t73 * t306) * t238;
t29 = -(t307 / 0.2e1 + t310 / 0.2e1) * t268 + t60;
t28 = -(t308 / 0.2e1 + t311 / 0.2e1) * t270 + t59;
t27 = -(t309 / 0.2e1 + t312 / 0.2e1) * t272 + t58;
t25 = -(t313 / 0.2e1 + t314 / 0.2e1) * t274 + t57;
t24 = -t159 * t56 + t64 * t80;
t23 = -t159 * t55 + t63 * t79;
t22 = -t159 * t54 + t62 * t78;
t21 = -t159 * t53 + t61 * t77;
t20 = -t56 * t80 + t64 * t76;
t19 = -t55 * t79 + t63 * t75;
t18 = -t54 * t78 + t62 * t74;
t17 = -t53 * t77 + t61 * t73;
t16 = (-pkin(3) * t316 + (pkin(3) * t32 + t60 * t243) * t60) * t327;
t15 = (-pkin(3) * t317 + (pkin(3) * t31 + t59 * t244) * t59) * t326;
t14 = (-pkin(3) * t318 + (pkin(3) * t30 + t58 * t245) * t58) * t325;
t13 = (-pkin(3) * t319 + (pkin(3) * t26 + t57 * t246) * t57) * t324;
t12 = ((-t255 * t29 * t276 - t235 * t32 - t237 * t60) * t236 * t60 + (pkin(3) + t243) * t316) * t327;
t11 = ((-t256 * t28 * t276 - t235 * t31 - t237 * t59) * t236 * t59 + (pkin(3) + t244) * t317) * t326;
t10 = ((-t257 * t27 * t276 - t235 * t30 - t237 * t58) * t236 * t58 + (pkin(3) + t245) * t318) * t325;
t9 = ((-t258 * t25 * t276 - t235 * t26 - t237 * t57) * t236 * t57 + (pkin(3) + t246) * t319) * t324;
t8 = t159 * t12 + t80 * t16 + (-t259 * pkin(2) * t60 ^ 2 + t263) * m(3);
t7 = t159 * t11 + t79 * t15 + (-t260 * pkin(2) * t59 ^ 2 + t264) * m(3);
t6 = t159 * t10 + t78 * t14 + (-t261 * pkin(2) * t58 ^ 2 + t265) * m(3);
t5 = t77 * t13 + t159 * t9 + (-t262 * pkin(2) * t57 ^ 2 + t266) * m(3);
t4 = t80 * t12 + t76 * t16 + ((-rSges(2,1) * t154 + rSges(2,2) * t158) * t189 + (rSges(2,1) * t158 + rSges(2,2) * t154) * t186) * m(2) + ((-0.2e1 * t239 * t259 * t29 - t154 * t189 + t158 * t186) * pkin(2) + t263) * m(3);
t3 = t79 * t11 + t75 * t15 + ((-rSges(2,1) * t153 + rSges(2,2) * t157) * t188 + (rSges(2,1) * t157 + rSges(2,2) * t153) * t185) * m(2) + ((-0.2e1 * t240 * t260 * t28 - t153 * t188 + t157 * t185) * pkin(2) + t264) * m(3);
t2 = t78 * t10 + t74 * t14 + ((-rSges(2,1) * t152 + rSges(2,2) * t156) * t187 + (rSges(2,1) * t156 + rSges(2,2) * t152) * t184) * m(2) + ((-0.2e1 * t241 * t261 * t27 - t152 * t187 + t156 * t184) * pkin(2) + t265) * m(3);
t1 = t73 * t13 + t77 * t9 + ((-rSges(2,1) * t151 + rSges(2,2) * t155) * t182 + (rSges(2,1) * t155 + rSges(2,2) * t151) * t181) * m(2) + ((-0.2e1 * t242 * t25 * t262 - t151 * t182 + t155 * t181) * pkin(2) + t266) * m(3);
t41 = [(-t135 * t215 - t210 * t136 - g(1) + t218) * m(4) + ((t123 * t40 * t92 + (t40 * t96 + t4) * t124) * t320 + (t121 * t38 * t91 + (t38 * t95 + t3) * t122) * t321 + (t119 * t36 * t90 + (t36 * t94 + t2) * t120) * t322 + (t117 * t34 * t89 + (t34 * t93 + t1) * t118) * t323 + (-(t8 * t88 + (-t88 * t269 + t80 * t294) * t251) * t320 - (t7 * t87 + (-t87 * t271 + t79 * t297) * t252) * t321 - (t6 * t86 + (-t86 * t273 + t78 * t300) * t253) * t322 - (t5 * t82 + (-t82 * t275 + t77 * t305) * t254) * t323) * t236) * t238; (-t210 * t135 + t136 * t215 - g(2) + t217) * m(4) + ((t124 * t39 * t96 + (t39 * t92 + t4) * t123) * t320 + (t122 * t37 * t95 + (t37 * t91 + t3) * t121) * t321 + (t120 * t35 * t94 + (t35 * t90 + t2) * t119) * t322 + (t118 * t33 * t93 + (t33 * t89 + t1) * t117) * t323 + (-(t8 * t85 + (-t85 * t269 + t80 * t295) * t251) * t320 - (t7 * t84 + (-t84 * t271 + t79 * t298) * t252) * t321 - (t6 * t83 + (-t83 * t273 + t78 * t301) * t253) * t322 - (t5 * t81 + (-t81 * t275 + t77 * t306) * t254) * t323) * t236) * t238; (-g(3) + xDDP(3)) * ((4 * m(1)) + 0.4e1 * m(2) + 0.4e1 * m(3) + m(4)); Icges(4,3) * t215 + t61 * t1 + t62 * t2 + t63 * t3 + t64 * t4 - t53 * t5 - t54 * t6 - t55 * t7 - t56 * t8 + (-t135 * t218 + t136 * t217 + (t223 ^ 2 + t224 ^ 2) * t215 + (g(1) * t223 - g(2) * t224) * t204 + (g(1) * t224 + g(2) * t223) * t203) * m(4) + ((t20 * t294 - t24 * t277) * t96 + (t20 * t295 - t24 * t278) * t92 + (t19 * t297 - t23 * t279) * t95 + (t19 * t298 - t23 * t280) * t91 + (t18 * t300 - t22 * t281) * t94 + (t18 * t301 - t22 * t282) * t90 + (t17 * t305 - t21 * t283) * t93 + (t17 * t306 - t21 * t284) * t89) * t238;];
tauX  = t41;
