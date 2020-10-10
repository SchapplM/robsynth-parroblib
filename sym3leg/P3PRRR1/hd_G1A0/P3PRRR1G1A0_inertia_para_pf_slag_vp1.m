% Calculate inertia matrix for parallel robot
% P3PRRR1G1A0
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRR1G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:47
% EndTime: 2020-03-09 21:14:47
% DurationCPUTime: 0.32s
% Computational Cost: add. (1633->86), mult. (1119->165), div. (180->5), fcn. (936->18), ass. (0->89)
t271 = pkin(7) + qJ(2,3);
t262 = qJ(3,3) + t271;
t248 = sin(t262);
t251 = cos(t262);
t256 = sin(t271);
t259 = cos(t271);
t308 = 0.1e1 / (t248 * t259 - t256 * t251);
t272 = pkin(7) + qJ(2,2);
t263 = qJ(3,2) + t272;
t249 = sin(t263);
t252 = cos(t263);
t257 = sin(t272);
t260 = cos(t272);
t307 = 0.1e1 / (t249 * t260 - t257 * t252);
t273 = pkin(7) + qJ(2,1);
t264 = qJ(3,1) + t273;
t250 = sin(t264);
t253 = cos(t264);
t258 = sin(t273);
t261 = cos(t273);
t306 = 0.1e1 / (t250 * t261 - t258 * t253);
t305 = m(3) * pkin(2);
t274 = legFrame(3,3);
t265 = sin(t274);
t268 = cos(t274);
t234 = t268 * t248 + t265 * t251;
t219 = pkin(2) * (t256 * t268 + t265 * t259) + t234 * pkin(3);
t304 = t219 * t308;
t275 = legFrame(2,3);
t266 = sin(t275);
t269 = cos(t275);
t236 = t269 * t249 + t266 * t252;
t220 = pkin(2) * (t257 * t269 + t266 * t260) + t236 * pkin(3);
t303 = t220 * t307;
t276 = legFrame(1,3);
t267 = sin(t276);
t270 = cos(t276);
t238 = t270 * t250 + t267 * t253;
t221 = pkin(2) * (t258 * t270 + t267 * t261) + t238 * pkin(3);
t302 = t221 * t306;
t235 = -t248 * t265 + t268 * t251;
t222 = -pkin(2) * (t256 * t265 - t268 * t259) + t235 * pkin(3);
t301 = t222 * t308;
t237 = -t249 * t266 + t269 * t252;
t223 = -pkin(2) * (t257 * t266 - t269 * t260) + t237 * pkin(3);
t300 = t223 * t307;
t239 = -t250 * t267 + t270 * t253;
t224 = -pkin(2) * (t258 * t267 - t270 * t261) + t239 * pkin(3);
t299 = t224 * t306;
t298 = t308 * t234;
t297 = t308 * t235;
t296 = t307 * t236;
t295 = t307 * t237;
t294 = t306 * t238;
t293 = t306 * t239;
t291 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t246 = t291 * m(3) + Icges(3,3);
t279 = 0.1e1 / pkin(3);
t292 = t246 * t279;
t283 = ((t256 * rSges(3,1) - rSges(3,2) * t259) * t248 + (rSges(3,1) * t259 + t256 * rSges(3,2)) * t251) * t305;
t216 = t283 + t246;
t290 = t216 * t308 * t279;
t282 = ((t257 * rSges(3,1) - rSges(3,2) * t260) * t249 + (rSges(3,1) * t260 + t257 * rSges(3,2)) * t252) * t305;
t217 = t282 + t246;
t289 = t217 * t307 * t279;
t281 = ((t258 * rSges(3,1) - rSges(3,2) * t261) * t250 + (rSges(3,1) * t261 + t258 * rSges(3,2)) * t253) * t305;
t218 = t281 + t246;
t288 = t218 * t306 * t279;
t287 = t308 * t292;
t286 = t307 * t292;
t285 = t306 * t292;
t284 = Icges(2,3) + Icges(3,3) + (pkin(2) ^ 2 + t291) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t280 = 0.1e1 / pkin(2);
t215 = 0.2e1 * t281 + t284;
t214 = 0.2e1 * t282 + t284;
t213 = 0.2e1 * t283 + t284;
t212 = (t218 * t293 - t224 * t285) * t280;
t211 = (t218 * t294 - t221 * t285) * t280;
t210 = (t217 * t295 - t223 * t286) * t280;
t209 = (t217 * t296 - t220 * t286) * t280;
t208 = (t216 * t297 - t222 * t287) * t280;
t207 = (t216 * t298 - t219 * t287) * t280;
t206 = (t215 * t293 - t224 * t288) * t280;
t205 = (t215 * t294 - t221 * t288) * t280;
t204 = (t214 * t295 - t223 * t289) * t280;
t203 = (t214 * t296 - t220 * t289) * t280;
t202 = (t213 * t297 - t222 * t290) * t280;
t201 = (t213 * t298 - t219 * t290) * t280;
t1 = [m(4) + (t202 * t297 + t204 * t295 + t206 * t293 + (-t208 * t301 - t210 * t300 - t212 * t299) * t279) * t280, (t202 * t298 + t204 * t296 + t206 * t294 + (-t208 * t304 - t210 * t303 - t212 * t302) * t279) * t280, 0; (t201 * t297 + t203 * t295 + t205 * t293 + (-t207 * t301 - t209 * t300 - t211 * t299) * t279) * t280, m(4) + (t201 * t298 + t203 * t296 + t205 * t294 + (-t207 * t304 - t209 * t303 - t211 * t302) * t279) * t280, 0; 0, 0, (3 * m(1)) + 0.3e1 * m(2) + 0.3e1 * m(3) + m(4);];
MX  = t1;
