% Calculate vector of inverse dynamics forces for parallel robot
% P4PRRRR1G3P1A0
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
% tauX [4x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:00:37
% EndTime: 2020-03-02 19:00:48
% DurationCPUTime: 10.65s
% Computational Cost: add. (5133->493), mult. (10727->947), div. (4160->18), fcn. (8222->34), ass. (0->329)
t171 = cos(qJ(3,4));
t148 = 0.1e1 / t171;
t187 = cos(qJ(3,3));
t156 = 0.1e1 / t187;
t189 = cos(qJ(3,2));
t160 = 0.1e1 / t189;
t191 = cos(qJ(3,1));
t164 = 0.1e1 / t191;
t218 = 0.1e1 / pkin(2);
t353 = t218 ^ 2;
t135 = (m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4));
t352 = 2 * t135;
t196 = m(2) * rSges(2,1);
t351 = rSges(3,3) * m(3);
t216 = rSges(3,2) ^ 2;
t217 = rSges(3,1) ^ 2;
t125 = (-t216 + t217) * m(3) - Icges(3,1) + Icges(3,2);
t350 = t125 / 0.2e1;
t169 = sin(qJ(3,4));
t333 = rSges(3,2) * t171;
t115 = rSges(3,1) * t169 + t333;
t349 = m(3) * t115;
t181 = sin(qJ(3,3));
t332 = rSges(3,2) * t187;
t116 = rSges(3,1) * t181 + t332;
t348 = m(3) * t116;
t183 = sin(qJ(3,2));
t331 = rSges(3,2) * t189;
t117 = rSges(3,1) * t183 + t331;
t347 = m(3) * t117;
t185 = sin(qJ(3,1));
t330 = rSges(3,2) * t191;
t118 = rSges(3,1) * t185 + t330;
t346 = m(3) * t118;
t197 = xDP(4);
t167 = t197 ^ 2;
t345 = m(4) * t167;
t177 = xDDP(4);
t344 = m(4) * t177;
t179 = xDDP(2);
t343 = m(4) * t179;
t180 = xDDP(1);
t342 = m(4) * t180;
t170 = sin(qJ(2,4));
t146 = 0.1e1 / t170;
t173 = legFrame(4,2);
t136 = sin(t173);
t140 = cos(t173);
t107 = g(1) * t136 + g(2) * t140;
t198 = xDP(3);
t282 = t198 * t218;
t252 = 0.2e1 * m(3) * t282;
t131 = rSges(3,2) * t252;
t132 = m(2) * rSges(2,2) - t351;
t151 = m(1) + m(2) + m(3);
t237 = rSges(3,1) * t252;
t147 = t171 ^ 2;
t150 = t148 / t147;
t168 = t198 ^ 2;
t291 = t168 / pkin(2) ^ 2;
t265 = t150 * t291;
t311 = t148 * t169;
t172 = cos(qJ(2,4));
t201 = xP(4);
t144 = sin(t201);
t145 = cos(t201);
t208 = koppelP(4,2);
t212 = koppelP(4,1);
t103 = -t144 * t208 + t145 * t212;
t199 = xDP(2);
t200 = xDP(1);
t289 = t169 * t198;
t149 = 0.1e1 / t171 ^ 2;
t309 = t149 * t172;
t99 = t144 * t212 + t145 * t208;
t42 = (-t289 * t309 + (-t140 * (-t197 * t99 + t200) + t136 * (t103 * t197 + t199)) * t148) * t218 * t146;
t321 = t172 * t42;
t292 = t168 * t218;
t41 = t42 ^ 2;
t33 = (-pkin(2) * t171 * t41 - t150 * t292) * t146;
t337 = rSges(3,1) * t171;
t37 = t149 * t291 + t41;
t236 = -rSges(3,2) * t169 + t337;
t81 = (t236 * m(3) + t196) * t172 - t170 * t132;
t267 = t148 * t282;
t290 = t169 * t171;
t308 = t149 * t218;
t9 = ((pkin(2) * t147 * t321 - t170 * t289) * t42 * t308 + (-t170 * t42 * t290 + t172 * t267) * t150 * t282) * t146;
t5 = -t81 * t9 - (t132 * t42 + t237 * t311 + t131) * t321 + (-t41 * t196 + (-t37 * t337 + (rSges(3,2) * t37 - t115 * t265) * t169) * m(3)) * t170 + (-t33 - t107) * t151;
t341 = t146 * t5;
t182 = sin(qJ(2,3));
t152 = 0.1e1 / t182;
t174 = legFrame(3,2);
t137 = sin(t174);
t141 = cos(t174);
t108 = g(1) * t137 + g(2) * t141;
t155 = t187 ^ 2;
t158 = t156 / t155;
t188 = cos(qJ(2,3));
t261 = t156 * t282;
t287 = t181 * t198;
t288 = t181 * t187;
t157 = 0.1e1 / t187 ^ 2;
t301 = t157 * t218;
t209 = koppelP(3,2);
t213 = koppelP(3,1);
t100 = t144 * t213 + t145 * t209;
t104 = -t144 * t209 + t145 * t213;
t302 = t157 * t188;
t46 = (-t287 * t302 + (-t141 * (-t100 * t197 + t200) + t137 * (t104 * t197 + t199)) * t156) * t218 * t152;
t320 = t188 * t46;
t12 = ((pkin(2) * t155 * t320 - t182 * t287) * t46 * t301 + (-t182 * t46 * t288 + t188 * t261) * t158 * t282) * t152;
t259 = t158 * t291;
t304 = t156 * t181;
t336 = rSges(3,1) * t187;
t43 = t46 ^ 2;
t34 = (-pkin(2) * t187 * t43 - t158 * t292) * t152;
t38 = t157 * t291 + t43;
t235 = -rSges(3,2) * t181 + t336;
t82 = (t235 * m(3) + t196) * t188 - t182 * t132;
t6 = -t82 * t12 - (t132 * t46 + t237 * t304 + t131) * t320 + (-t43 * t196 + (-t38 * t336 + (rSges(3,2) * t38 - t116 * t259) * t181) * m(3)) * t182 + (-t34 - t108) * t151;
t340 = t152 * t6;
t184 = sin(qJ(2,2));
t153 = 0.1e1 / t184;
t159 = t189 ^ 2;
t162 = t160 / t159;
t190 = cos(qJ(2,2));
t258 = t160 * t282;
t285 = t183 * t198;
t286 = t183 * t189;
t161 = 0.1e1 / t189 ^ 2;
t297 = t161 * t218;
t210 = koppelP(2,2);
t214 = koppelP(2,1);
t101 = t144 * t214 + t145 * t210;
t105 = -t144 * t210 + t145 * t214;
t175 = legFrame(2,2);
t138 = sin(t175);
t142 = cos(t175);
t298 = t161 * t190;
t47 = (-t285 * t298 + (-t142 * (-t101 * t197 + t200) + t138 * (t105 * t197 + t199)) * t160) * t218 * t153;
t319 = t190 * t47;
t10 = ((pkin(2) * t159 * t319 - t184 * t285) * t47 * t297 + (-t184 * t47 * t286 + t190 * t258) * t162 * t282) * t153;
t109 = g(1) * t138 + g(2) * t142;
t256 = t162 * t291;
t300 = t160 * t183;
t335 = rSges(3,1) * t189;
t44 = t47 ^ 2;
t35 = (-pkin(2) * t189 * t44 - t162 * t292) * t153;
t39 = t161 * t291 + t44;
t234 = -rSges(3,2) * t183 + t335;
t83 = (t234 * m(3) + t196) * t190 - t184 * t132;
t7 = -t83 * t10 - (t132 * t47 + t237 * t300 + t131) * t319 + (-t44 * t196 + (-t39 * t335 + (rSges(3,2) * t39 - t117 * t256) * t183) * m(3)) * t184 + (-t35 - t109) * t151;
t339 = t153 * t7;
t186 = sin(qJ(2,1));
t154 = 0.1e1 / t186;
t163 = t191 ^ 2;
t166 = t164 / t163;
t192 = cos(qJ(2,1));
t255 = t164 * t282;
t283 = t185 * t198;
t284 = t185 * t191;
t165 = 0.1e1 / t191 ^ 2;
t293 = t165 * t218;
t211 = koppelP(1,2);
t215 = koppelP(1,1);
t102 = t144 * t215 + t145 * t211;
t106 = -t144 * t211 + t145 * t215;
t176 = legFrame(1,2);
t139 = sin(t176);
t143 = cos(t176);
t294 = t165 * t192;
t48 = (-t283 * t294 + (-t143 * (-t102 * t197 + t200) + t139 * (t106 * t197 + t199)) * t164) * t218 * t154;
t318 = t192 * t48;
t11 = ((pkin(2) * t163 * t318 - t186 * t283) * t48 * t293 + (-t186 * t48 * t284 + t192 * t255) * t166 * t282) * t154;
t110 = g(1) * t139 + g(2) * t143;
t253 = t166 * t291;
t296 = t164 * t185;
t334 = rSges(3,1) * t191;
t45 = t48 ^ 2;
t36 = (-pkin(2) * t191 * t45 - t166 * t292) * t154;
t40 = t165 * t291 + t45;
t233 = -rSges(3,2) * t185 + t334;
t84 = (t233 * m(3) + t196) * t192 - t186 * t132;
t8 = -t84 * t11 - (t132 * t48 + t237 * t296 + t131) * t318 + (-t45 * t196 + (-t40 * t334 + (rSges(3,2) * t40 - t118 * t253) * t185) * m(3)) * t186 + (-t36 - t110) * t151;
t338 = t154 * t8;
t329 = t146 * (t103 * t177 - t167 * t99 + t179);
t328 = t146 * (-t103 * t167 - t177 * t99 + t180);
t327 = t152 * (-t100 * t167 + t104 * t177 + t179);
t326 = t152 * (-t100 * t177 - t104 * t167 + t180);
t325 = t153 * (-t101 * t167 + t105 * t177 + t179);
t324 = t153 * (-t101 * t177 - t105 * t167 + t180);
t323 = t154 * (-t102 * t167 + t106 * t177 + t179);
t322 = t154 * (-t102 * t177 - t106 * t167 + t180);
t317 = t115 * t170;
t316 = t116 * t182;
t315 = t117 * t184;
t314 = t118 * t186;
t281 = t216 + t217;
t130 = t281 * m(3) + Icges(3,3);
t313 = t130 * t218;
t312 = t146 * t169;
t310 = t148 * t218;
t307 = t152 * t181;
t306 = t153 * t183;
t305 = t154 * t185;
t303 = t156 * t218;
t299 = t160 * t218;
t295 = t164 * t218;
t280 = m(3) * t317;
t279 = m(3) * t316;
t278 = m(3) * t315;
t277 = m(3) * t314;
t276 = t136 * t310;
t275 = t137 * t303;
t274 = t138 * t299;
t273 = t139 * t295;
t272 = t140 * t310;
t271 = t141 * t303;
t270 = t142 * t299;
t269 = t143 * t295;
t268 = t146 * t310;
t266 = t172 * t308;
t264 = t152 * t303;
t263 = t153 * t299;
t262 = t154 * t295;
t260 = t188 * t301;
t257 = t190 * t297;
t254 = t192 * t293;
t251 = t136 * t268;
t250 = t137 * t264;
t249 = t138 * t263;
t248 = t139 * t262;
t247 = t140 * t268;
t246 = t141 * t264;
t245 = t142 * t263;
t244 = t143 * t262;
t243 = t169 * t265;
t242 = t181 * t259;
t241 = t183 * t256;
t240 = t185 * t253;
t133 = rSges(3,2) * t351 - Icges(3,6);
t239 = -t133 * t282 / 0.4e1;
t134 = rSges(3,1) * t351 - Icges(3,5);
t238 = t134 * t282 / 0.4e1;
t232 = t266 * t312;
t231 = t260 * t307;
t230 = t257 * t306;
t229 = t254 * t305;
t111 = g(1) * t140 - g(2) * t136;
t228 = t107 * t170 + t111 * t172;
t112 = g(1) * t141 - g(2) * t137;
t227 = t108 * t182 + t112 * t188;
t113 = g(1) * t142 - g(2) * t138;
t226 = t109 * t184 + t113 * t190;
t114 = g(1) * t143 - g(2) * t139;
t225 = t110 * t186 + t114 * t192;
t224 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + ((2 * rSges(3,3) ^ 2) + t281) * m(3) / 0.2e1;
t207 = 0.2e1 * qJ(3,1);
t206 = 0.2e1 * qJ(3,2);
t205 = 0.2e1 * qJ(3,3);
t204 = 0.2e1 * qJ(3,4);
t203 = rSges(4,1);
t202 = rSges(4,2);
t195 = rSges(3,2) * g(3);
t178 = xDDP(3);
t98 = -t144 * t202 + t145 * t203;
t97 = t144 * t203 + t145 * t202;
t96 = t139 * t186 + t143 * t192;
t95 = -t139 * t192 + t143 * t186;
t94 = t138 * t184 + t142 * t190;
t93 = -t138 * t190 + t142 * t184;
t92 = t137 * t182 + t141 * t188;
t91 = -t137 * t188 + t141 * t182;
t90 = t136 * t170 + t140 * t172;
t89 = -t136 * t172 + t140 * t170;
t88 = -t133 * t191 - t134 * t185;
t87 = -t133 * t189 - t134 * t183;
t86 = -t133 * t187 - t134 * t181;
t85 = -t133 * t171 - t134 * t169;
t72 = cos(t207) * t350 - t135 * sin(t207) + t224;
t71 = cos(t206) * t350 - t135 * sin(t206) + t224;
t70 = cos(t205) * t350 - t135 * sin(t205) + t224;
t69 = cos(t204) * t350 - t135 * sin(t204) + t224;
t68 = (t102 * t143 + t106 * t139) * t262;
t67 = (t101 * t142 + t105 * t138) * t263;
t66 = (t100 * t141 + t104 * t137) * t264;
t65 = (t103 * t136 + t140 * t99) * t268;
t64 = (t151 * t96 - t84 * t269) * t154;
t63 = (t151 * t95 + t84 * t273) * t154;
t62 = (t151 * t94 - t83 * t270) * t153;
t61 = (t151 * t93 + t83 * t274) * t153;
t60 = (t151 * t92 - t82 * t271) * t152;
t59 = (t151 * t91 + t82 * t275) * t152;
t58 = (t151 * t90 - t81 * t272) * t146;
t57 = (t151 * t89 + t81 * t276) * t146;
t56 = (-t102 * t96 + t106 * t95) * t154;
t55 = (-t101 * t94 + t105 * t93) * t153;
t54 = (-t100 * t92 + t104 * t91) * t152;
t53 = (t103 * t89 - t90 * t99) * t146;
t52 = -t277 * t295 + (t151 * t164 - t84 * t254) * t305;
t51 = -t278 * t299 + (t151 * t160 - t83 * t257) * t306;
t50 = -t279 * t303 + (t151 * t156 - t82 * t260) * t307;
t49 = -t280 * t310 + (t148 * t151 - t81 * t266) * t312;
t32 = (-t72 * t269 + t84 * t96) * t154;
t31 = (t72 * t273 + t84 * t95) * t154;
t30 = (-t71 * t270 + t83 * t94) * t153;
t29 = (t71 * t274 + t83 * t93) * t153;
t28 = (-t70 * t271 + t82 * t92) * t152;
t27 = (t70 * t275 + t82 * t91) * t152;
t26 = (-t69 * t272 + t81 * t90) * t146;
t25 = (t69 * t276 + t81 * t89) * t146;
t24 = t88 * t295 + (t164 * t84 - t72 * t254) * t305;
t23 = t87 * t299 + (t160 * t83 - t71 * t257) * t306;
t22 = t86 * t303 + (t156 * t82 - t70 * t260) * t307;
t21 = t85 * t310 + (t148 * t81 - t69 * t266) * t312;
t20 = t151 * t56 + t68 * t84;
t19 = t151 * t55 + t67 * t83;
t18 = t151 * t54 + t66 * t82;
t17 = t151 * t53 + t65 * t81;
t16 = t56 * t84 + t68 * t72;
t15 = t55 * t83 + t67 * t71;
t14 = t54 * t82 + t66 * t70;
t13 = t53 * t81 + t65 * t69;
t4 = -t84 * t36 - t72 * t11 + t88 * t240 - 0.4e1 * ((t48 * t185 * t350 + t164 * t238) * t191 + t239 * t296 + (t163 - 0.1e1 / 0.2e1) * t48 * t135) * t255 + (m(2) * (-rSges(2,1) * t110 + rSges(2,2) * t114) + (-t114 * rSges(3,3) - t233 * t110) * m(3)) * t192 + t186 * (m(2) * (rSges(2,1) * t114 + rSges(2,2) * t110) + (-t110 * rSges(3,3) + t233 * t114) * m(3));
t3 = -t83 * t35 - t71 * t10 + t87 * t241 - 0.4e1 * ((t47 * t183 * t350 + t160 * t238) * t189 + t239 * t300 + (t159 - 0.1e1 / 0.2e1) * t47 * t135) * t258 + (m(2) * (-rSges(2,1) * t109 + rSges(2,2) * t113) + (-t113 * rSges(3,3) - t234 * t109) * m(3)) * t190 + t184 * (m(2) * (rSges(2,1) * t113 + rSges(2,2) * t109) + (-t109 * rSges(3,3) + t234 * t113) * m(3));
t2 = -t82 * t34 - t70 * t12 + t86 * t242 - 0.4e1 * ((t46 * t181 * t350 + t156 * t238) * t187 + t239 * t304 + (t155 - 0.1e1 / 0.2e1) * t46 * t135) * t261 + (m(2) * (-rSges(2,1) * t108 + rSges(2,2) * t112) + (-t112 * rSges(3,3) - t235 * t108) * m(3)) * t188 + t182 * (m(2) * (rSges(2,1) * t112 + rSges(2,2) * t108) + (-t108 * rSges(3,3) + t235 * t112) * m(3));
t1 = -t81 * t33 - t69 * t9 + t85 * t243 - 0.4e1 * ((t42 * t169 * t350 + t148 * t238) * t171 + t239 * t311 + (t147 - 0.1e1 / 0.2e1) * t42 * t135) * t267 + (m(2) * (-rSges(2,1) * t107 + rSges(2,2) * t111) + (-t111 * rSges(3,3) - t236 * t107) * m(3)) * t172 + t170 * (m(2) * (rSges(2,1) * t111 + rSges(2,2) * t107) + (-t107 * rSges(3,3) + t236 * t111) * m(3));
t73 = [(-t32 * t269 + t64 * t96) * t322 + (t32 * t273 + t64 * t95) * t323 + t96 * t338 - t4 * t244 + (-t30 * t270 + t62 * t94) * t324 + (t30 * t274 + t62 * t93) * t325 + t94 * t339 - t3 * t245 + (-t28 * t271 + t60 * t92) * t326 + (t28 * t275 + t60 * t91) * t327 + t92 * t340 - t2 * t246 + (-t26 * t272 + t58 * t90) * t328 + (t26 * t276 + t58 * t89) * t329 + t90 * t341 - t1 * t247 + t342 - t97 * t344 - t98 * t345 - m(4) * g(1) + ((-t88 * t244 - t96 * t346) * t295 + (t164 * t64 - t32 * t254) * t305 + (-t87 * t245 - t94 * t347) * t299 + (t160 * t62 - t30 * t257) * t306 + (-t86 * t246 - t92 * t348) * t303 + (t156 * t60 - t28 * t260) * t307 + (-t85 * t247 - t90 * t349) * t310 + (t148 * t58 - t26 * t266) * t312) * t178; (-t31 * t269 + t63 * t96) * t322 + (t31 * t273 + t63 * t95) * t323 + t95 * t338 + t4 * t248 + (-t29 * t270 + t61 * t94) * t324 + (t29 * t274 + t61 * t93) * t325 + t93 * t339 + t3 * t249 + (-t27 * t271 + t59 * t92) * t326 + (t27 * t275 + t59 * t91) * t327 + t91 * t340 + t2 * t250 + (-t25 * t272 + t57 * t90) * t328 + (t25 * t276 + t57 * t89) * t329 + t89 * t341 + t1 * t251 + t343 + t98 * t344 - t97 * t345 - m(4) * g(2) + ((t88 * t248 - t95 * t346) * t295 + (t164 * t63 - t31 * t254) * t305 + (t87 * t249 - t93 * t347) * t299 + (t160 * t61 - t29 * t257) * t306 + (t86 * t250 - t91 * t348) * t303 + (t156 * t59 - t27 * t260) * t307 + (t85 * t251 - t89 * t349) * t310 + (t148 * t57 - t25 * t266) * t312) * t178; -t1 * t232 - t4 * t229 - t3 * t230 - t2 * t231 + (-t88 * t11 + t130 * t240 + (t125 * t284 + t163 * t352 + Icges(3,4)) * t45 + (t36 * t314 + t185 * t195 + t225 * t330 + (-rSges(3,2) * t45 - g(3) * t191 + t185 * t225) * rSges(3,1)) * m(3)) * t295 + (-t87 * t10 + t130 * t241 + (t125 * t286 + t159 * t352 + Icges(3,4)) * t44 + (t35 * t315 + t183 * t195 + t226 * t331 + (-rSges(3,2) * t44 - g(3) * t189 + t183 * t226) * rSges(3,1)) * m(3)) * t299 + (-t86 * t12 + t130 * t242 + (t125 * t288 + t155 * t352 + Icges(3,4)) * t43 + (t34 * t316 + t181 * t195 + t227 * t332 + (-rSges(3,2) * t43 - g(3) * t187 + t181 * t227) * rSges(3,1)) * m(3)) * t303 + (-t85 * t9 + t130 * t243 + (t125 * t290 + t147 * t352 + Icges(3,4)) * t41 + (t33 * t317 + t169 * t195 + t228 * t333 + (-rSges(3,2) * t41 - g(3) * t171 + t169 * t228) * rSges(3,1)) * m(3)) * t310 - m(4) * g(3) + (-t24 * t269 + t52 * t96) * t322 + (t24 * t273 + t52 * t95) * t323 + (-t23 * t270 + t51 * t94) * t324 + (t23 * t274 + t51 * t93) * t325 + (-t22 * t271 + t50 * t92) * t326 + (t22 * t275 + t50 * t91) * t327 + (-t21 * t272 + t49 * t90) * t328 + (t21 * t276 + t49 * t89) * t329 + t296 * t338 + t300 * t339 + t304 * t340 + t311 * t341 + (-t21 * t232 + ((-t169 * t349 + t313) * t310 + (-t353 * t85 * t309 + t49) * t312) * t148 - t24 * t229 + ((-t185 * t346 + t313) * t295 + (-t353 * t88 * t294 + t52) * t305) * t164 - t23 * t230 + ((-t183 * t347 + t313) * t299 + (-t353 * t87 * t298 + t51) * t306) * t160 - t22 * t231 + ((-t181 * t348 + t313) * t303 + (-t353 * t86 * t302 + t50) * t307) * t156 + m(4)) * t178; (-t16 * t269 + t20 * t96) * t322 + (t16 * t273 + t20 * t95) * t323 + t56 * t8 + t68 * t4 + (-t15 * t270 + t19 * t94) * t324 + (t15 * t274 + t19 * t93) * t325 + t55 * t7 + t67 * t3 + (-t14 * t271 + t18 * t92) * t326 + (t14 * t275 + t18 * t91) * t327 + t54 * t6 + t66 * t2 + (-t13 * t272 + t17 * t90) * t328 + (t13 * t276 + t17 * t89) * t329 + t53 * t5 + t65 * t1 - t97 * t342 + t98 * t343 + (Icges(4,3) + m(4) * (t202 ^ 2 + t203 ^ 2)) * t177 + m(4) * ((g(1) * t203 + g(2) * t202) * t144 + (g(1) * t202 - g(2) * t203) * t145) + ((-t56 * t277 + t68 * t88) * t295 + (-t16 * t254 + t164 * t20) * t305 + (-t55 * t278 + t67 * t87) * t299 + (-t15 * t257 + t160 * t19) * t306 + (-t54 * t279 + t66 * t86) * t303 + (-t14 * t260 + t156 * t18) * t307 + (-t53 * t280 + t65 * t85) * t310 + (-t13 * t266 + t148 * t17) * t312) * t178;];
tauX  = t73;
