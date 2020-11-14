% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
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
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:58:59
% EndTime: 2020-08-06 19:59:04
% DurationCPUTime: 4.76s
% Computational Cost: add. (17061->456), mult. (28683->725), div. (3201->12), fcn. (21807->44), ass. (0->313)
t193 = pkin(4) + qJ(3,3);
t179 = 0.1e1 / t193;
t220 = xDP(3);
t221 = xDP(2);
t222 = xDP(1);
t196 = legFrame(3,2);
t158 = sin(t196);
t161 = cos(t196);
t203 = sin(qJ(1,3));
t209 = cos(qJ(1,3));
t189 = cos(pkin(5));
t137 = t189 * pkin(2) + pkin(1);
t208 = cos(qJ(2,3));
t188 = sin(pkin(5));
t202 = sin(qJ(2,3));
t285 = t188 * t202;
t236 = pkin(2) * t285 - t137 * t208;
t83 = -t209 * t193 - t236 * t203;
t331 = pkin(2) * t188;
t95 = t202 * t137 + t208 * t331;
t55 = t95 * t158 + t83 * t161;
t56 = -t83 * t158 + t95 * t161;
t175 = qJ(2,3) + pkin(5);
t155 = cos(t175);
t168 = t208 * pkin(1);
t355 = pkin(2) * t155 + t168;
t86 = t203 * t193 + t209 * t355;
t40 = (t220 * t86 + t221 * t56 + t222 * t55) * t179;
t368 = 0.2e1 * t40;
t194 = pkin(4) + qJ(3,2);
t180 = 0.1e1 / t194;
t197 = legFrame(2,2);
t159 = sin(t197);
t162 = cos(t197);
t205 = sin(qJ(1,2));
t211 = cos(qJ(1,2));
t210 = cos(qJ(2,2));
t204 = sin(qJ(2,2));
t284 = t188 * t204;
t235 = pkin(2) * t284 - t137 * t210;
t84 = -t211 * t194 - t235 * t205;
t96 = t204 * t137 + t210 * t331;
t57 = t96 * t159 + t84 * t162;
t58 = -t84 * t159 + t96 * t162;
t176 = qJ(2,2) + pkin(5);
t156 = cos(t176);
t169 = t210 * pkin(1);
t354 = pkin(2) * t156 + t169;
t87 = t205 * t194 + t211 * t354;
t41 = (t220 * t87 + t221 * t58 + t222 * t57) * t180;
t367 = 0.2e1 * t41;
t195 = pkin(4) + qJ(3,1);
t181 = 0.1e1 / t195;
t198 = legFrame(1,2);
t160 = sin(t198);
t163 = cos(t198);
t207 = sin(qJ(1,1));
t213 = cos(qJ(1,1));
t212 = cos(qJ(2,1));
t206 = sin(qJ(2,1));
t283 = t188 * t206;
t234 = pkin(2) * t283 - t137 * t212;
t85 = -t213 * t195 - t234 * t207;
t97 = t206 * t137 + t212 * t331;
t59 = t97 * t160 + t85 * t163;
t60 = -t85 * t160 + t97 * t163;
t177 = qJ(2,1) + pkin(5);
t157 = cos(t177);
t170 = t212 * pkin(1);
t353 = pkin(2) * t157 + t170;
t88 = t207 * t195 + t213 * t353;
t42 = (t220 * t88 + t221 * t60 + t222 * t59) * t181;
t366 = 0.2e1 * t42;
t276 = 0.2e1 * pkin(1);
t275 = t189 * t276;
t231 = pkin(2) ^ 2;
t232 = pkin(1) ^ 2;
t277 = -t231 - t232;
t116 = pkin(2) * t275 - t277;
t106 = 0.1e1 / t355;
t73 = (t158 * t222 + t161 * t221) * t106;
t70 = t73 ^ 2;
t364 = t70 * t116;
t107 = 0.1e1 / t354;
t74 = (t159 * t222 + t162 * t221) * t107;
t71 = t74 ^ 2;
t363 = t71 * t116;
t108 = 0.1e1 / t353;
t75 = (t160 * t222 + t163 * t221) * t108;
t72 = t75 ^ 2;
t362 = t72 * t116;
t149 = sin(t175);
t165 = t202 * pkin(1);
t361 = t149 * rSges(3,1) + t155 * rSges(3,2) + t165;
t360 = pkin(2) * t149 + t165;
t150 = sin(t176);
t166 = t204 * pkin(1);
t359 = t150 * rSges(3,1) + t156 * rSges(3,2) + t166;
t358 = pkin(2) * t150 + t166;
t151 = sin(t177);
t167 = t206 * pkin(1);
t357 = t151 * rSges(3,1) + t157 * rSges(3,2) + t167;
t356 = pkin(2) * t151 + t167;
t346 = m(3) * rSges(3,2);
t273 = t188 * t346;
t98 = -m(2) * rSges(2,1) + t273 + (-rSges(3,1) * t189 - pkin(1)) * m(3);
t348 = m(2) * rSges(2,2);
t99 = t348 + (rSges(3,1) * t188 + rSges(3,2) * t189) * m(3);
t352 = -t99 * t206 - t98 * t212;
t351 = -t99 * t204 - t98 * t210;
t350 = -t99 * t202 - t98 * t208;
t349 = m(1) * rSges(1,1);
t347 = m(2) * rSges(2,3);
t89 = t155 * rSges(3,1) - t149 * rSges(3,2) + t168;
t345 = m(3) * t89;
t90 = t156 * rSges(3,1) - t150 * rSges(3,2) + t169;
t344 = m(3) * t90;
t91 = t157 * rSges(3,1) - t151 * rSges(3,2) + t170;
t343 = m(3) * t91;
t228 = rSges(2,2) ^ 2;
t230 = rSges(2,1) ^ 2;
t109 = m(3) * t232 + (-t228 + t230) * m(2) + Icges(2,2) - Icges(2,1);
t342 = t109 / 0.2e1;
t227 = rSges(3,2) ^ 2;
t229 = rSges(3,1) ^ 2;
t118 = (-t227 + t229) * m(3) - Icges(3,1) + Icges(3,2);
t341 = t118 / 0.2e1;
t340 = m(3) * t179;
t339 = m(3) * t180;
t338 = m(3) * t181;
t190 = qJ(3,3) + rSges(3,3);
t330 = t190 * m(3);
t191 = qJ(3,2) + rSges(3,3);
t329 = t191 * m(3);
t192 = qJ(3,1) + rSges(3,3);
t328 = t192 * m(3);
t266 = t158 * t331;
t269 = t161 * t331;
t288 = t158 * t137;
t297 = t137 * t161;
t61 = (-t203 * t288 + t269) * t208 + t202 * (t203 * t266 + t297);
t92 = 0.1e1 / t236;
t327 = t61 * t92;
t265 = t159 * t331;
t268 = t162 * t331;
t287 = t159 * t137;
t296 = t137 * t162;
t62 = (-t205 * t287 + t268) * t210 + t204 * (t205 * t265 + t296);
t93 = 0.1e1 / t235;
t326 = t62 * t93;
t264 = t160 * t331;
t267 = t163 * t331;
t286 = t160 * t137;
t295 = t137 * t163;
t63 = (-t207 * t286 + t267) * t212 + t206 * (t207 * t264 + t295);
t94 = 0.1e1 / t234;
t325 = t63 * t94;
t64 = (t203 * t297 + t266) * t208 + (-t203 * t269 + t288) * t202;
t324 = t64 * t92;
t65 = (t205 * t296 + t265) * t210 + (-t205 * t268 + t287) * t204;
t323 = t65 * t93;
t66 = (t207 * t295 + t264) * t212 + (-t207 * t267 + t286) * t206;
t322 = t66 * t94;
t321 = t89 * t92;
t320 = t90 * t93;
t319 = t91 * t94;
t31 = (t209 * t220 - (t221 * t61 + t222 * t64) * t92) * t179;
t318 = t31 / 0.2e1;
t32 = (t211 * t220 - (t221 * t62 + t222 * t65) * t93) * t180;
t317 = t32 / 0.2e1;
t33 = (t213 * t220 - (t221 * t63 + t222 * t66) * t94) * t181;
t316 = t33 / 0.2e1;
t309 = t106 * t360;
t308 = t106 * t158;
t307 = t106 * t161;
t306 = t107 * t358;
t305 = t107 * t159;
t304 = t107 * t162;
t303 = t108 * t356;
t302 = t108 * t160;
t301 = t108 * t163;
t224 = 0.2e1 * qJ(2,3);
t182 = sin(t224);
t300 = t109 * t182;
t225 = 0.2e1 * qJ(2,2);
t183 = sin(t225);
t299 = t109 * t183;
t226 = 0.2e1 * qJ(2,1);
t184 = sin(t226);
t298 = t109 * t184;
t138 = 0.2e1 * t175;
t124 = cos(t138);
t142 = -rSges(3,1) * t346 + Icges(3,4);
t294 = t142 * t124;
t139 = 0.2e1 * t176;
t125 = cos(t139);
t293 = t142 * t125;
t140 = 0.2e1 * t177;
t126 = cos(t140);
t292 = t142 * t126;
t144 = rSges(2,1) * t348 - Icges(2,4);
t185 = cos(t224);
t291 = t144 * t185;
t186 = cos(t225);
t290 = t144 * t186;
t187 = cos(t226);
t289 = t144 * t187;
t282 = -rSges(2,1) * t347 + Icges(2,5);
t172 = pkin(5) + t224;
t152 = cos(t172);
t281 = t152 + t189;
t173 = pkin(5) + t225;
t153 = cos(t173);
t280 = t153 + t189;
t174 = pkin(5) + t226;
t154 = cos(t174);
t279 = t154 + t189;
t278 = t228 + t230;
t127 = rSges(3,2) * t330 - Icges(3,6);
t130 = rSges(3,1) * t330 - Icges(3,5);
t143 = rSges(2,2) * t347 - Icges(2,6);
t239 = pkin(1) * t330 - t282;
t46 = (t127 * t188 - t130 * t189 - t239) * t202 - t208 * (t127 * t189 + t130 * t188 + t143);
t272 = t179 * t46 * t92;
t128 = rSges(3,2) * t329 - Icges(3,6);
t131 = rSges(3,1) * t329 - Icges(3,5);
t238 = pkin(1) * t329 - t282;
t47 = (t128 * t188 - t131 * t189 - t238) * t204 - t210 * (t128 * t189 + t131 * t188 + t143);
t271 = t180 * t47 * t93;
t129 = rSges(3,2) * t328 - Icges(3,6);
t132 = rSges(3,1) * t328 - Icges(3,5);
t237 = pkin(1) * t328 - t282;
t48 = (t129 * t188 - t132 * t189 - t237) * t206 - t212 * (t129 * t189 + t132 * t188 + t143);
t270 = t181 * t48 * t94;
t263 = t46 * t308;
t262 = t47 * t305;
t261 = t48 * t302;
t260 = t46 * t307;
t259 = t47 * t304;
t258 = t48 * t301;
t257 = m(3) * t276;
t256 = pkin(2) * t276;
t249 = t227 / 0.2e1 + t229 / 0.2e1 + t232 / 0.2e1;
t146 = sin(t172);
t248 = rSges(3,1) * t146 + rSges(3,2) * t152;
t147 = sin(t173);
t247 = rSges(3,1) * t147 + rSges(3,2) * t153;
t148 = sin(t174);
t246 = rSges(3,1) * t148 + rSges(3,2) * t154;
t103 = t161 * g(1) - t158 * g(2);
t242 = g(3) * t209 + t103 * t203;
t104 = t162 * g(1) - t159 * g(2);
t241 = g(3) * t211 + t104 * t205;
t105 = t163 * g(1) - t160 * g(2);
t240 = g(3) * t213 + t105 * t207;
t136 = pkin(1) * t273;
t233 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,2) / 0.2e1 - t136 + (0.2e1 * rSges(2,3) ^ 2 + t278) * m(2) / 0.2e1;
t201 = xDDP(1);
t200 = xDDP(2);
t199 = xDDP(3);
t178 = g(3) * t349;
t141 = m(1) * rSges(1,2) - t347;
t133 = t141 * g(3);
t123 = sin(t140);
t122 = sin(t139);
t121 = sin(t138);
t102 = t160 * g(1) + t163 * g(2);
t101 = t159 * g(1) + t162 * g(2);
t100 = t158 * g(1) + t161 * g(2);
t79 = -0.2e1 * t136 + t278 * m(2) + Icges(2,3) + Icges(3,3) + (rSges(3,1) * t275 + t227 + t229 + t232) * m(3);
t69 = (-t213 * t91 + t88) * t338;
t68 = (-t211 * t90 + t87) * t339;
t67 = (-t209 * t89 + t86) * t340;
t45 = t126 * t341 + t187 * t342 + (t192 ^ 2 + (t279 * rSges(3,1) - rSges(3,2) * t148) * pkin(1) + t249) * m(3) + t142 * t123 - t144 * t184 + t233;
t44 = t125 * t341 + t186 * t342 + (t191 ^ 2 + (t280 * rSges(3,1) - rSges(3,2) * t147) * pkin(1) + t249) * m(3) + t142 * t122 - t144 * t183 + t233;
t43 = t124 * t341 + t185 * t342 + (t190 ^ 2 + (t281 * rSges(3,1) - rSges(3,2) * t146) * pkin(1) + t249) * m(3) + t142 * t121 - t144 * t182 + t233;
t39 = (t66 * t319 + t59) * t338;
t38 = (t65 * t320 + t57) * t339;
t37 = (t64 * t321 + t55) * t340;
t36 = (t63 * t319 + t60) * t338;
t35 = (t62 * t320 + t58) * t339;
t34 = (t61 * t321 + t56) * t340;
t30 = (t213 * t45 - t88 * t343) * t181;
t29 = (t211 * t44 - t87 * t344) * t180;
t28 = (t209 * t43 - t86 * t345) * t179;
t27 = -t66 * t270 + t79 * t302;
t26 = -t65 * t271 + t79 * t305;
t25 = -t64 * t272 + t79 * t308;
t24 = -t63 * t270 + t79 * t301;
t23 = -t62 * t271 + t79 * t304;
t22 = -t61 * t272 + t79 * t307;
t21 = t261 + (-t45 * t322 - t59 * t343) * t181;
t20 = t262 + (-t44 * t323 - t57 * t344) * t180;
t19 = t263 + (-t43 * t324 - t55 * t345) * t179;
t18 = t258 + (-t45 * t325 - t60 * t343) * t181;
t17 = t259 + (-t44 * t326 - t58 * t344) * t180;
t16 = t260 + (-t43 * t327 - t56 * t345) * t179;
t15 = (-0.1e1 / (t170 + (t189 * t212 - t283) * pkin(2)) * t362 + (t234 * t33 + t366) * t33) * t181;
t14 = (-0.1e1 / (t169 + (t189 * t210 - t284) * pkin(2)) * t363 + (t235 * t32 + t367) * t32) * t180;
t13 = (-0.1e1 / (t168 + (t189 * t208 - t285) * pkin(2)) * t364 + (t236 * t31 + t368) * t31) * t179;
t12 = (-t362 + ((-t231 * t126 - t187 * t232 - t279 * t256 + t277) * t316 + t353 * t366 + 0.2e1 * (-t195 * t316 + t356 * t75) * t195) * t33) * t181;
t11 = (-t363 + ((-t231 * t125 - t186 * t232 - t280 * t256 + t277) * t317 + t354 * t367 + 0.2e1 * (-t194 * t317 + t358 * t74) * t194) * t32) * t180;
t10 = (-t364 + ((-t231 * t124 - t185 * t232 - t281 * t256 + t277) * t318 + t355 * t368 + 0.2e1 * (-t193 * t318 + t360 * t73) * t193) * t31) * t179;
t9 = -t48 * t15 + t79 * t72 * t303 + (t102 * t99 - t240 * t98) * t206 + (t102 * t98 + t240 * t99) * t212 + ((t123 * t341 - t292 + t298 / 0.2e1 + t289) * t33 + (t246 * t33 * pkin(1) - 0.2e1 * t357 * t42) * m(3)) * t33;
t8 = -t47 * t14 + t79 * t71 * t306 + (t101 * t99 - t241 * t98) * t204 + (t101 * t98 + t241 * t99) * t210 + ((t122 * t341 - t293 + t299 / 0.2e1 + t290) * t32 + (t247 * t32 * pkin(1) - 0.2e1 * t359 * t41) * m(3)) * t32;
t7 = -t46 * t13 + t79 * t70 * t309 + (t100 * t99 - t242 * t98) * t202 + (t100 * t98 + t242 * t99) * t208 + ((t121 * t341 - t294 + t300 / 0.2e1 + t291) * t31 + (t248 * t31 * pkin(1) - 0.2e1 * t361 * t40) * m(3)) * t31;
t6 = (-g(3) * t207 + t105 * t213 + t91 * t15 - t12 + (-t33 * t192 + 0.2e1 * t357 * t75) * t33) * m(3);
t5 = (-g(3) * t205 + t104 * t211 + t90 * t14 - t11 + (-t32 * t191 + 0.2e1 * t359 * t74) * t32) * m(3);
t4 = (-g(3) * t203 + t103 * t209 + t89 * t13 - t10 + (-t31 * t190 + 0.2e1 * t361 * t73) * t31) * m(3);
t3 = -t45 * t15 + t12 * t343 + (-g(3) * t328 + t133) * t213 + (g(3) * t352 + t178) * t207 + ((-t349 - t352) * t213 + (t141 - t328) * t207) * t105 + (t129 * t151 - t132 * t157 + t143 * t206 - t237 * t212 + t48 * t303) * t72 + (t328 * t366 + (-t118 * t123 - t246 * t257 - 0.2e1 * t289 + 0.2e1 * t292 - t298) * t75) * t33;
t2 = -t44 * t14 + t11 * t344 + (-g(3) * t329 + t133) * t211 + (g(3) * t351 + t178) * t205 + ((-t349 - t351) * t211 + (t141 - t329) * t205) * t104 + (t128 * t150 - t131 * t156 + t143 * t204 - t238 * t210 + t47 * t306) * t71 + (t329 * t367 + (-t118 * t122 - t247 * t257 - 0.2e1 * t290 + 0.2e1 * t293 - t299) * t74) * t32;
t1 = -t43 * t13 + t10 * t345 + (-g(3) * t330 + t133) * t209 + (g(3) * t350 + t178) * t203 + ((-t349 - t350) * t209 + (t141 - t330) * t203) * t103 + (t127 * t149 - t130 * t155 + t143 * t202 - t239 * t208 + t46 * t309) * t70 + (t330 * t368 + (-t118 * t121 - t248 * t257 - 0.2e1 * t291 + 0.2e1 * t294 - t300) * t73) * t31;
t49 = [t7 * t308 + t8 * t305 + t9 * t302 - m(4) * g(1) + (t25 * t307 + t26 * t304 + t27 * t301) * t200 + (t25 * t308 + t26 * t305 + t27 * t302 + m(4)) * t201 + ((-t21 * t322 + t39 * t59) * t201 + (-t21 * t325 + t39 * t60) * t200 + (t21 * t213 + t39 * t88) * t199 - t3 * t322 + t59 * t6) * t181 + ((-t20 * t323 + t38 * t57) * t201 + (-t20 * t326 + t38 * t58) * t200 + (t20 * t211 + t38 * t87) * t199 - t2 * t323 + t57 * t5) * t180 + ((-t19 * t324 + t37 * t55) * t201 + (-t19 * t327 + t37 * t56) * t200 + (t19 * t209 + t37 * t86) * t199 - t1 * t324 + t55 * t4) * t179; t7 * t307 + t8 * t304 + t9 * t301 - m(4) * g(2) + (t22 * t308 + t23 * t305 + t24 * t302) * t201 + (t22 * t307 + t23 * t304 + t24 * t301 + m(4)) * t200 + ((-t18 * t322 + t36 * t59) * t201 + (-t18 * t325 + t36 * t60) * t200 + (t18 * t213 + t36 * t88) * t199 - t3 * t325 + t60 * t6) * t181 + ((-t17 * t323 + t35 * t57) * t201 + (-t17 * t326 + t35 * t58) * t200 + (t17 * t211 + t35 * t87) * t199 - t2 * t326 + t58 * t5) * t180 + ((-t16 * t324 + t34 * t55) * t201 + (-t16 * t327 + t34 * t56) * t200 + (t16 * t209 + t34 * t86) * t199 - t1 * t327 + t56 * t4) * t179; (-g(3) + t199) * m(4) + ((t213 * t261 - t30 * t322 + t59 * t69) * t201 + (t213 * t258 - t30 * t325 + t60 * t69) * t200 + (t213 * t30 + t69 * t88) * t199 + t213 * t3 + t88 * t6) * t181 + ((t211 * t262 - t29 * t323 + t57 * t68) * t201 + (t211 * t259 - t29 * t326 + t58 * t68) * t200 + (t211 * t29 + t68 * t87) * t199 + t211 * t2 + t87 * t5) * t180 + ((t209 * t263 - t28 * t324 + t55 * t67) * t201 + (t209 * t260 - t28 * t327 + t56 * t67) * t200 + (t209 * t28 + t67 * t86) * t199 + t209 * t1 + t86 * t4) * t179;];
tauX  = t49;
