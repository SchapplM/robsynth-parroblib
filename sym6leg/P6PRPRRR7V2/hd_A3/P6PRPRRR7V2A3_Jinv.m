% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d2,d4,theta1,theta3]';
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [6x6]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-17 04:59
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6PRPRRR7V2A3_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(10,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6PRPRRR7V2A3_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6PRPRRR7V2A3_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PRPRRR7V2A3_Jinv: pkin has to be [10x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6PRPRRR7V2A3_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PRPRRR7V2A3_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-17 04:52:52
% EndTime: 2019-05-17 04:52:57
% DurationCPUTime: 6.20s
% Computational Cost: add. (1962->346), mult. (4224->682), div. (36->12), fcn. (4236->64), ass. (0->325)
t213 = cos(pkin(5));
t279 = sin(pkin(4)) * t213;
t210 = cos(pkin(10));
t352 = sin(pkin(10));
t359 = sin(pkin(6)) * pkin(8);
t141 = pkin(3) * t210 + t352 * t359 + pkin(2);
t241 = sin(qJ(2,1));
t244 = cos(qJ(2,1));
t140 = (-pkin(3) * t352 + t210 * t359) * t213;
t161 = cos(pkin(6)) * pkin(8);
t160 = t161 + qJ(3,1);
t208 = sin(pkin(5));
t272 = t208 * t160 + t140;
t360 = t241 * t141 - t272 * t244;
t371 = t360 * t279;
t240 = sin(qJ(2,2));
t243 = cos(qJ(2,2));
t159 = t161 + qJ(3,2);
t273 = t208 * t159 + t140;
t361 = t240 * t141 - t273 * t243;
t370 = t361 * t279;
t239 = sin(qJ(2,3));
t242 = cos(qJ(2,3));
t158 = t161 + qJ(3,3);
t274 = t208 * t158 + t140;
t362 = t239 * t141 - t274 * t242;
t369 = t362 * t279;
t229 = sin(qJ(2,4));
t232 = cos(qJ(2,4));
t157 = t161 + qJ(3,4);
t275 = t208 * t157 + t140;
t363 = t229 * t141 - t275 * t232;
t368 = t363 * t279;
t228 = sin(qJ(2,5));
t231 = cos(qJ(2,5));
t156 = t161 + qJ(3,5);
t276 = t208 * t156 + t140;
t364 = t228 * t141 - t276 * t231;
t367 = t364 * t279;
t227 = sin(qJ(2,6));
t230 = cos(qJ(2,6));
t155 = t161 + qJ(3,6);
t277 = t208 * t155 + t140;
t365 = t227 * t141 - t277 * t230;
t366 = t365 * t279;
t214 = cos(pkin(4));
t358 = t214 * (t141 * t230 + t277 * t227);
t357 = t214 * (t141 * t231 + t276 * t228);
t356 = t214 * (t141 * t232 + t275 * t229);
t355 = t214 * (t141 * t242 + t274 * t239);
t354 = t214 * (t141 * t243 + t273 * t240);
t353 = t214 * (t141 * t244 + t272 * t241);
t246 = xP(5);
t200 = sin(t246);
t203 = cos(t246);
t248 = koppelP(6,3);
t245 = xP(6);
t199 = sin(t245);
t202 = cos(t245);
t254 = koppelP(6,2);
t260 = koppelP(6,1);
t271 = -t199 * t254 + t202 * t260;
t109 = t200 * t248 + t271 * t203;
t206 = sin(pkin(9));
t351 = t109 * t206;
t211 = cos(pkin(9));
t350 = t109 * t211;
t249 = koppelP(5,3);
t255 = koppelP(5,2);
t261 = koppelP(5,1);
t270 = -t199 * t255 + t202 * t261;
t110 = t200 * t249 + t270 * t203;
t349 = t110 * t206;
t348 = t110 * t211;
t250 = koppelP(4,3);
t256 = koppelP(4,2);
t262 = koppelP(4,1);
t269 = -t199 * t256 + t202 * t262;
t111 = t200 * t250 + t269 * t203;
t347 = t111 * t206;
t346 = t111 * t211;
t251 = koppelP(3,3);
t257 = koppelP(3,2);
t263 = koppelP(3,1);
t268 = -t199 * t257 + t202 * t263;
t112 = t200 * t251 + t268 * t203;
t345 = t112 * t206;
t344 = t112 * t211;
t252 = koppelP(2,3);
t258 = koppelP(2,2);
t264 = koppelP(2,1);
t267 = -t199 * t258 + t202 * t264;
t113 = t200 * t252 + t203 * t267;
t343 = t113 * t206;
t342 = t113 * t211;
t253 = koppelP(1,3);
t259 = koppelP(1,2);
t265 = koppelP(1,1);
t266 = -t199 * t259 + t202 * t265;
t114 = t200 * t253 + t203 * t266;
t341 = t114 * t206;
t340 = t114 * t211;
t221 = legFrame(6,1);
t169 = sin(t221);
t233 = legFrame(6,2);
t187 = sin(t233);
t339 = t169 * t187;
t222 = legFrame(5,1);
t170 = sin(t222);
t234 = legFrame(5,2);
t188 = sin(t234);
t338 = t170 * t188;
t223 = legFrame(4,1);
t171 = sin(t223);
t235 = legFrame(4,2);
t189 = sin(t235);
t337 = t171 * t189;
t224 = legFrame(3,1);
t172 = sin(t224);
t236 = legFrame(3,2);
t190 = sin(t236);
t336 = t172 * t190;
t225 = legFrame(2,1);
t173 = sin(t225);
t237 = legFrame(2,2);
t191 = sin(t237);
t335 = t173 * t191;
t226 = legFrame(1,1);
t174 = sin(t226);
t238 = legFrame(1,2);
t192 = sin(t238);
t334 = t174 * t192;
t181 = cos(t221);
t333 = t181 * t187;
t182 = cos(t222);
t332 = t182 * t188;
t183 = cos(t223);
t331 = t183 * t189;
t184 = cos(t224);
t330 = t184 * t190;
t185 = cos(t225);
t329 = t185 * t191;
t186 = cos(t226);
t328 = t186 * t192;
t327 = t187 * t206;
t326 = t187 * t211;
t325 = t188 * t206;
t324 = t188 * t211;
t323 = t189 * t206;
t322 = t189 * t211;
t321 = t190 * t206;
t320 = t190 * t211;
t319 = t191 * t206;
t318 = t191 * t211;
t317 = t192 * t206;
t316 = t192 * t211;
t247 = xP(4);
t201 = sin(t247);
t315 = t201 * t254;
t314 = t201 * t255;
t313 = t201 * t256;
t312 = t201 * t257;
t311 = t201 * t258;
t310 = t201 * t259;
t309 = t201 * t260;
t308 = t201 * t261;
t307 = t201 * t262;
t306 = t201 * t263;
t305 = t201 * t264;
t304 = t201 * t265;
t303 = t203 * t248;
t302 = t203 * t249;
t301 = t203 * t250;
t300 = t203 * t251;
t299 = t203 * t252;
t298 = t203 * t253;
t204 = cos(t247);
t297 = t204 * t254;
t296 = t204 * t255;
t295 = t204 * t256;
t294 = t204 * t257;
t293 = t204 * t258;
t292 = t204 * t259;
t291 = t204 * t260;
t290 = t204 * t261;
t289 = t204 * t262;
t288 = t204 * t263;
t287 = t204 * t264;
t286 = t204 * t265;
t115 = t208 * t140;
t278 = -t115 - t161;
t220 = legFrame(1,3);
t219 = legFrame(2,3);
t218 = legFrame(3,3);
t217 = legFrame(4,3);
t216 = legFrame(5,3);
t215 = legFrame(6,3);
t205 = t213 ^ 2;
t198 = cos(t238);
t197 = cos(t237);
t196 = cos(t236);
t195 = cos(t235);
t194 = cos(t234);
t193 = cos(t233);
t180 = cos(t220);
t179 = cos(t219);
t178 = cos(t218);
t177 = cos(t217);
t176 = cos(t216);
t175 = cos(t215);
t168 = sin(t220);
t167 = sin(t219);
t166 = sin(t218);
t165 = sin(t217);
t164 = sin(t216);
t163 = sin(t215);
t154 = t160 * t205;
t153 = t159 * t205;
t152 = t158 * t205;
t151 = t157 * t205;
t150 = t156 * t205;
t149 = t155 * t205;
t147 = t199 * t265 + t202 * t259;
t146 = t199 * t264 + t202 * t258;
t145 = t199 * t263 + t202 * t257;
t144 = t199 * t262 + t202 * t256;
t143 = t199 * t261 + t202 * t255;
t142 = t199 * t260 + t202 * t254;
t139 = -t206 * t168 + t211 * t180;
t138 = -t206 * t167 + t211 * t179;
t137 = -t206 * t166 + t211 * t178;
t136 = -t206 * t165 + t211 * t177;
t135 = -t206 * t164 + t211 * t176;
t134 = -t206 * t163 + t211 * t175;
t133 = t211 * t168 + t206 * t180;
t132 = t211 * t167 + t206 * t179;
t131 = t211 * t166 + t206 * t178;
t130 = t211 * t165 + t206 * t177;
t129 = t211 * t164 + t206 * t176;
t128 = t211 * t163 + t206 * t175;
t102 = t200 * t266 - t298;
t101 = t200 * t267 - t299;
t100 = t268 * t200 - t300;
t99 = t269 * t200 - t301;
t98 = t270 * t200 - t302;
t97 = t271 * t200 - t303;
t78 = t102 * t204 - t201 * t147;
t77 = t101 * t204 - t201 * t146;
t76 = t100 * t204 - t201 * t145;
t75 = -t201 * t144 + t99 * t204;
t74 = -t201 * t143 + t98 * t204;
t73 = -t201 * t142 + t97 * t204;
t72 = t201 * t102 + t204 * t147;
t71 = t201 * t101 + t204 * t146;
t70 = t201 * t100 + t204 * t145;
t69 = t204 * t144 + t201 * t99;
t68 = t204 * t143 + t201 * t98;
t67 = t204 * t142 + t201 * t97;
t66 = (-t204 * t298 + (t200 * t286 - t310) * t202 - t199 * (t200 * t292 + t304)) * t198 - t114 * t328;
t65 = (-t201 * t298 + (t200 * t304 + t292) * t202 + t199 * (-t200 * t310 + t286)) * t198 - t114 * t334;
t64 = (-t204 * t299 + (t200 * t287 - t311) * t202 - t199 * (t200 * t293 + t305)) * t197 - t113 * t329;
t63 = (-t201 * t299 + (t200 * t305 + t293) * t202 + t199 * (-t200 * t311 + t287)) * t197 - t113 * t335;
t62 = (-t204 * t300 + (t200 * t288 - t312) * t202 - t199 * (t200 * t294 + t306)) * t196 - t112 * t330;
t61 = (-t201 * t300 + (t200 * t306 + t294) * t202 + t199 * (-t200 * t312 + t288)) * t196 - t112 * t336;
t60 = (-t204 * t301 + (t200 * t289 - t313) * t202 - t199 * (t200 * t295 + t307)) * t195 - t111 * t331;
t59 = (-t201 * t301 + (t200 * t307 + t295) * t202 + t199 * (-t200 * t313 + t289)) * t195 - t111 * t337;
t58 = (-t204 * t302 + (t200 * t290 - t314) * t202 - t199 * (t200 * t296 + t308)) * t194 - t110 * t332;
t57 = (-t201 * t302 + (t200 * t308 + t296) * t202 + t199 * (-t200 * t314 + t290)) * t194 - t110 * t338;
t56 = (-t204 * t303 + (t200 * t291 - t315) * t202 - t199 * (t200 * t297 + t309)) * t193 - t109 * t333;
t55 = (-t201 * t303 + (t200 * t309 + t297) * t202 + t199 * (-t200 * t315 + t291)) * t193 - t109 * t339;
t54 = 0.1e1 / ((t154 - t115 - t160) * t214 + t371);
t53 = 0.1e1 / ((t153 - t115 - t159) * t214 + t370);
t52 = 0.1e1 / ((t152 - t115 - t158) * t214 + t369);
t51 = 0.1e1 / ((t151 - t115 - t157) * t214 + t368);
t50 = 0.1e1 / ((t150 - t115 - t156) * t214 + t367);
t49 = 0.1e1 / ((t149 - t115 - t155) * t214 + t366);
t48 = -t206 * t78 + t72 * t316;
t47 = t211 * t78 + t72 * t317;
t46 = -t206 * t77 + t71 * t318;
t45 = t211 * t77 + t71 * t319;
t44 = -t206 * t76 + t70 * t320;
t43 = t211 * t76 + t70 * t321;
t42 = -t206 * t75 + t69 * t322;
t41 = t211 * t75 + t69 * t323;
t40 = -t206 * t74 + t68 * t324;
t39 = t211 * t74 + t68 * t325;
t38 = -t206 * t73 + t67 * t326;
t37 = t211 * t73 + t67 * t327;
t36 = t72 * t206 + t78 * t316;
t35 = t71 * t206 + t77 * t318;
t34 = t70 * t206 + t76 * t320;
t33 = t69 * t206 + t75 * t322;
t32 = t68 * t206 + t74 * t324;
t31 = t67 * t206 + t73 * t326;
t30 = -t72 * t211 + t78 * t317;
t29 = -t71 * t211 + t77 * t319;
t28 = -t70 * t211 + t76 * t321;
t27 = -t69 * t211 + t75 * t323;
t26 = -t68 * t211 + t74 * t325;
t25 = -t67 * t211 + t73 * t327;
t24 = t174 * t341 + t66 * t211;
t23 = -t186 * t341 + t65 * t211;
t22 = t173 * t343 + t64 * t211;
t21 = -t185 * t343 + t63 * t211;
t20 = t172 * t345 + t62 * t211;
t19 = -t184 * t345 + t61 * t211;
t18 = t171 * t347 + t60 * t211;
t17 = -t183 * t347 + t59 * t211;
t16 = t170 * t349 + t58 * t211;
t15 = -t182 * t349 + t57 * t211;
t14 = t169 * t351 + t56 * t211;
t13 = -t181 * t351 + t55 * t211;
t12 = t186 * t340 + t65 * t206;
t11 = t185 * t342 + t63 * t206;
t10 = t184 * t344 + t61 * t206;
t9 = t183 * t346 + t59 * t206;
t8 = t182 * t348 + t57 * t206;
t7 = t181 * t350 + t55 * t206;
t6 = t174 * t340 - t206 * t66;
t5 = t173 * t342 - t206 * t64;
t4 = t172 * t344 - t206 * t62;
t3 = t171 * t346 - t206 * t60;
t2 = t170 * t348 - t206 * t58;
t1 = t169 * t350 - t206 * t56;
t79 = [-t198 * (-t133 * t360 + t139 * t353) / ((-qJ(3,1) + t154 + t278) * t214 + t371), (-(t133 * t186 + t139 * t334) * t353 + (t133 * t334 - t139 * t186) * t360) * t54, ((-t174 * t133 + t139 * t328) * t353 - (t133 * t328 + t174 * t139) * t360) * t54, (-((t168 * t47 - t48 * t180) * t186 + (-t168 * t30 + t36 * t180) * t174) * t353 + t360 * ((-t48 * t168 - t47 * t180) * t186 + (t168 * t36 + t30 * t180) * t174)) * t54, ((t6 * t168 + t24 * t180) * t353 - t360 * (t24 * t168 - t6 * t180)) * t54, ((-t12 * t168 + t23 * t180) * t353 - (t12 * t180 + t23 * t168) * t360) * t54; -t197 * (-t132 * t361 + t138 * t354) / ((-qJ(3,2) + t153 + t278) * t214 + t370), (-(t132 * t185 + t138 * t335) * t354 + (t132 * t335 - t138 * t185) * t361) * t53, ((-t173 * t132 + t138 * t329) * t354 - (t132 * t329 + t173 * t138) * t361) * t53, (-((t167 * t45 - t46 * t179) * t185 + (-t167 * t29 + t35 * t179) * t173) * t354 + t361 * ((-t46 * t167 - t45 * t179) * t185 + (t167 * t35 + t29 * t179) * t173)) * t53, ((t5 * t167 + t22 * t179) * t354 - t361 * (t22 * t167 - t5 * t179)) * t53, ((-t11 * t167 + t21 * t179) * t354 - (t11 * t179 + t21 * t167) * t361) * t53; -t196 * (-t131 * t362 + t137 * t355) / ((-qJ(3,3) + t152 + t278) * t214 + t369), (-(t131 * t184 + t137 * t336) * t355 + (t131 * t336 - t137 * t184) * t362) * t52, ((-t172 * t131 + t137 * t330) * t355 - (t131 * t330 + t172 * t137) * t362) * t52, (-((t166 * t43 - t44 * t178) * t184 + (-t166 * t28 + t34 * t178) * t172) * t355 + t362 * ((-t44 * t166 - t43 * t178) * t184 + (t166 * t34 + t28 * t178) * t172)) * t52, ((t4 * t166 + t20 * t178) * t355 - t362 * (t20 * t166 - t4 * t178)) * t52, ((-t10 * t166 + t19 * t178) * t355 - (t10 * t178 + t19 * t166) * t362) * t52; -t195 * (-t130 * t363 + t136 * t356) / ((-qJ(3,4) + t151 + t278) * t214 + t368), (-(t130 * t183 + t136 * t337) * t356 + (t130 * t337 - t136 * t183) * t363) * t51, ((-t171 * t130 + t136 * t331) * t356 - (t130 * t331 + t171 * t136) * t363) * t51, (-((t165 * t41 - t42 * t177) * t183 + (-t165 * t27 + t33 * t177) * t171) * t356 + t363 * ((-t42 * t165 - t41 * t177) * t183 + (t165 * t33 + t27 * t177) * t171)) * t51, ((t3 * t165 + t18 * t177) * t356 - (t18 * t165 - t3 * t177) * t363) * t51, ((-t9 * t165 + t17 * t177) * t356 - (t17 * t165 + t9 * t177) * t363) * t51; -t194 * (-t129 * t364 + t135 * t357) / ((-qJ(3,5) + t150 + t278) * t214 + t367), (-(t129 * t182 + t135 * t338) * t357 + (t129 * t338 - t135 * t182) * t364) * t50, ((-t170 * t129 + t135 * t332) * t357 - (t129 * t332 + t170 * t135) * t364) * t50, (-((t164 * t39 - t40 * t176) * t182 + (-t164 * t26 + t32 * t176) * t170) * t357 + t364 * ((-t40 * t164 - t39 * t176) * t182 + (t164 * t32 + t26 * t176) * t170)) * t50, ((t16 * t176 + t2 * t164) * t357 - (t16 * t164 - t2 * t176) * t364) * t50, ((t15 * t176 - t8 * t164) * t357 - (t15 * t164 + t8 * t176) * t364) * t50; -t193 * (-t128 * t365 + t134 * t358) / ((-qJ(3,6) + t149 + t278) * t214 + t366), (-(t128 * t181 + t134 * t339) * t358 + (t128 * t339 - t134 * t181) * t365) * t49, ((-t169 * t128 + t134 * t333) * t358 - (t128 * t333 + t169 * t134) * t365) * t49, (-((t163 * t37 - t38 * t175) * t181 + (-t163 * t25 + t31 * t175) * t169) * t358 + t365 * ((-t38 * t163 - t37 * t175) * t181 + (t163 * t31 + t25 * t175) * t169)) * t49, ((t1 * t163 + t14 * t175) * t358 - (-t1 * t175 + t14 * t163) * t365) * t49, ((t13 * t175 - t7 * t163) * t358 - (t13 * t163 + t7 * t175) * t365) * t49;];
Jinv  = t79;
