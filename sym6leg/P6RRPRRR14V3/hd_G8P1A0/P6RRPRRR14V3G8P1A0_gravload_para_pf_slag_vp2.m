% Calculate Gravitation load for parallel robot
% P6RRPRRR14V3G8P1A0
% Use Code from Maple symbolic Code Generation
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [6x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-12 23:36
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P6RRPRRR14V3G8P1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(1,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G8P1A0_gravload_para_pf_slag_vp2: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G8P1A0_gravload_para_pf_slag_vp2: qJ has to be [3x6] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G8P1A0_gravload_para_pf_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6RRPRRR14V3G8P1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6RRPRRR14V3G8P1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P6RRPRRR14V3G8P1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G8P1A0_gravload_para_pf_slag_vp2: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G8P1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 1
% StartTime: 2020-03-12 23:33:09
% EndTime: 2020-03-12 23:33:10
% DurationCPUTime: 0.55s
% Computational Cost: add. (4451->475), mult. (9156->1036), div. (162->12), fcn. (9905->66), ass. (0->362)
unknown=NaN(6,1);
t1 = 0.1e1 / qJ(3,1);
t2 = sin(qJ(1,1));
t3 = cos(legFrame(1,3));
t5 = cos(qJ(1,1));
t6 = sin(legFrame(1,3));
t8 = t3 * t2 + t6 * t5;
t10 = cos(legFrame(1,2));
t11 = sin(qJ(2,1));
t12 = 0.1e1 / t11;
t16 = sin(legFrame(1,1));
t17 = sin(legFrame(1,2));
t18 = t17 * t16;
t20 = cos(legFrame(1,1));
t24 = t17 * t20;
t29 = -g(1) * t6 * t10 + g(2) * (-t6 * t18 + t3 * t20) + g(3) * (t3 * t16 + t6 * t24);
t31 = m(3) * qJ(3,1) - mrSges(2,2) + mrSges(3,3);
t32 = t31 * t29;
t34 = mrSges(2,1) + mrSges(3,1);
t35 = t34 * t29;
t36 = cos(qJ(2,1));
t38 = mrSges(1,2) - mrSges(3,2) - mrSges(2,3);
t49 = g(1) * t3 * t10 + g(2) * (t3 * t18 + t6 * t20) + g(3) * (t6 * t16 - t3 * t24);
t54 = t31 * t49;
t56 = t34 * t49;
t62 = t5 * (-mrSges(1,1) * t29 - t11 * t32 - t36 * t35 + t49 * t38) + (mrSges(1,1) * t49 + t11 * t54 + t38 * t29 + t36 * t56) * t2;
t67 = -t2 * t6 + t5 * t3;
t71 = t10 * t67 * t36 + t17 * t11;
t76 = t10 * t16;
t78 = t10 * t20;
t80 = g(1) * t17 - g(2) * t76 + g(3) * t78;
t89 = t36 * (-t2 * t32 - t34 * t80 - t5 * t54) - (-t2 * t35 + t31 * t80 - t5 * t56) * t11;
t94 = t10 * t67 * t11 - t36 * t17;
t101 = t29 * t2 * t11 + t49 * t11 * t5 - t80 * t36;
t103 = 0.1e1 / qJ(3,2);
t104 = sin(qJ(1,2));
t105 = cos(legFrame(2,3));
t107 = cos(qJ(1,2));
t108 = sin(legFrame(2,3));
t110 = t105 * t104 + t108 * t107;
t112 = cos(legFrame(2,2));
t113 = sin(qJ(2,2));
t114 = 0.1e1 / t113;
t118 = sin(legFrame(2,1));
t119 = sin(legFrame(2,2));
t120 = t119 * t118;
t122 = cos(legFrame(2,1));
t126 = t119 * t122;
t131 = -g(1) * t108 * t112 + g(2) * (t105 * t122 - t108 * t120) + g(3) * (t105 * t118 + t108 * t126);
t133 = m(3) * qJ(3,2) - mrSges(2,2) + mrSges(3,3);
t134 = t133 * t131;
t136 = t34 * t131;
t137 = cos(qJ(2,2));
t149 = g(1) * t105 * t112 + g(2) * (t105 * t120 + t108 * t122) + g(3) * (-t105 * t126 + t108 * t118);
t154 = t133 * t149;
t156 = t34 * t149;
t162 = t107 * (-mrSges(1,1) * t131 - t113 * t134 - t137 * t136 + t149 * t38) + (mrSges(1,1) * t149 + t113 * t154 + t38 * t131 + t137 * t156) * t104;
t167 = -t104 * t108 + t107 * t105;
t171 = t112 * t167 * t137 + t119 * t113;
t176 = t112 * t118;
t178 = t112 * t122;
t180 = g(1) * t119 - g(2) * t176 + g(3) * t178;
t189 = t137 * (-t104 * t134 - t107 * t154 - t34 * t180) - (-t104 * t136 - t107 * t156 + t133 * t180) * t113;
t194 = t112 * t167 * t113 - t137 * t119;
t201 = t131 * t104 * t113 + t149 * t113 * t107 - t180 * t137;
t203 = 0.1e1 / qJ(3,3);
t204 = sin(qJ(1,3));
t205 = cos(legFrame(3,3));
t207 = cos(qJ(1,3));
t208 = sin(legFrame(3,3));
t210 = t205 * t204 + t208 * t207;
t212 = cos(legFrame(3,2));
t213 = sin(qJ(2,3));
t214 = 0.1e1 / t213;
t218 = sin(legFrame(3,1));
t219 = sin(legFrame(3,2));
t220 = t219 * t218;
t222 = cos(legFrame(3,1));
t226 = t219 * t222;
t231 = -g(1) * t208 * t212 + g(2) * (t205 * t222 - t208 * t220) + g(3) * (t205 * t218 + t208 * t226);
t233 = m(3) * qJ(3,3) - mrSges(2,2) + mrSges(3,3);
t234 = t233 * t231;
t236 = t34 * t231;
t237 = cos(qJ(2,3));
t249 = g(1) * t205 * t212 + g(2) * (t205 * t220 + t208 * t222) + g(3) * (-t205 * t226 + t208 * t218);
t254 = t233 * t249;
t256 = t34 * t249;
t262 = t207 * (-mrSges(1,1) * t231 - t213 * t234 - t237 * t236 + t249 * t38) + (mrSges(1,1) * t249 + t213 * t254 + t38 * t231 + t237 * t256) * t204;
t267 = -t204 * t208 + t207 * t205;
t271 = t212 * t267 * t237 + t219 * t213;
t276 = t212 * t218;
t278 = t212 * t222;
t280 = g(1) * t219 - g(2) * t276 + g(3) * t278;
t289 = t237 * (-t204 * t234 - t207 * t254 - t34 * t280) - (-t204 * t236 - t207 * t256 + t233 * t280) * t213;
t294 = t212 * t267 * t213 - t237 * t219;
t301 = t231 * t204 * t213 + t249 * t213 * t207 - t280 * t237;
t303 = 0.1e1 / qJ(3,4);
t304 = sin(qJ(1,4));
t305 = cos(legFrame(4,3));
t307 = cos(qJ(1,4));
t308 = sin(legFrame(4,3));
t310 = t305 * t304 + t308 * t307;
t312 = cos(legFrame(4,2));
t313 = sin(qJ(2,4));
t314 = 0.1e1 / t313;
t318 = sin(legFrame(4,1));
t319 = sin(legFrame(4,2));
t320 = t319 * t318;
t322 = cos(legFrame(4,1));
t326 = t319 * t322;
t331 = -g(1) * t308 * t312 + g(2) * (t305 * t322 - t308 * t320) + g(3) * (t305 * t318 + t308 * t326);
t333 = m(3) * qJ(3,4) - mrSges(2,2) + mrSges(3,3);
t334 = t333 * t331;
t336 = t34 * t331;
t337 = cos(qJ(2,4));
t349 = g(1) * t305 * t312 + g(2) * (t305 * t320 + t308 * t322) + g(3) * (-t305 * t326 + t308 * t318);
t354 = t333 * t349;
t356 = t34 * t349;
t362 = t307 * (-mrSges(1,1) * t331 - t313 * t334 - t337 * t336 + t349 * t38) + (mrSges(1,1) * t349 + t313 * t354 + t38 * t331 + t337 * t356) * t304;
t367 = -t304 * t308 + t307 * t305;
t371 = t312 * t367 * t337 + t319 * t313;
t376 = t312 * t318;
t378 = t312 * t322;
t380 = g(1) * t319 - g(2) * t376 + g(3) * t378;
t389 = t337 * (-t304 * t334 - t307 * t354 - t34 * t380) - (-t304 * t336 - t307 * t356 + t333 * t380) * t313;
t394 = t312 * t367 * t313 - t337 * t319;
t401 = t331 * t304 * t313 + t349 * t313 * t307 - t380 * t337;
t403 = 0.1e1 / qJ(3,5);
t404 = sin(qJ(1,5));
t405 = cos(legFrame(5,3));
t407 = cos(qJ(1,5));
t408 = sin(legFrame(5,3));
t410 = t405 * t404 + t408 * t407;
t412 = cos(legFrame(5,2));
t413 = sin(qJ(2,5));
t414 = 0.1e1 / t413;
t418 = sin(legFrame(5,1));
t419 = sin(legFrame(5,2));
t420 = t419 * t418;
t422 = cos(legFrame(5,1));
t426 = t419 * t422;
t431 = -g(1) * t408 * t412 + g(2) * (t405 * t422 - t408 * t420) + g(3) * (t405 * t418 + t408 * t426);
t433 = m(3) * qJ(3,5) - mrSges(2,2) + mrSges(3,3);
t434 = t433 * t431;
t436 = t34 * t431;
t437 = cos(qJ(2,5));
t449 = g(1) * t405 * t412 + g(2) * (t405 * t420 + t408 * t422) + g(3) * (-t405 * t426 + t408 * t418);
t454 = t433 * t449;
t456 = t34 * t449;
t462 = t407 * (-mrSges(1,1) * t431 + t449 * t38 - t413 * t434 - t437 * t436) + (mrSges(1,1) * t449 + t38 * t431 + t413 * t454 + t437 * t456) * t404;
t467 = -t404 * t408 + t407 * t405;
t471 = t412 * t467 * t437 + t419 * t413;
t476 = t412 * t418;
t478 = t412 * t422;
t480 = g(1) * t419 - g(2) * t476 + g(3) * t478;
t489 = t437 * (-t34 * t480 - t404 * t434 - t407 * t454) - (-t404 * t436 - t407 * t456 + t433 * t480) * t413;
t494 = t412 * t467 * t413 - t437 * t419;
t501 = t431 * t404 * t413 + t449 * t413 * t407 - t480 * t437;
t503 = 0.1e1 / qJ(3,6);
t504 = sin(qJ(1,6));
t505 = cos(legFrame(6,3));
t507 = cos(qJ(1,6));
t508 = sin(legFrame(6,3));
t510 = t505 * t504 + t508 * t507;
t512 = cos(legFrame(6,2));
t513 = sin(qJ(2,6));
t514 = 0.1e1 / t513;
t518 = sin(legFrame(6,1));
t519 = sin(legFrame(6,2));
t520 = t519 * t518;
t522 = cos(legFrame(6,1));
t526 = t519 * t522;
t531 = -g(1) * t508 * t512 + g(2) * (t505 * t522 - t508 * t520) + g(3) * (t505 * t518 + t508 * t526);
t533 = m(3) * qJ(3,6) - mrSges(2,2) + mrSges(3,3);
t534 = t533 * t531;
t536 = t34 * t531;
t537 = cos(qJ(2,6));
t549 = g(1) * t505 * t512 + g(2) * (t505 * t520 + t508 * t522) + g(3) * (-t505 * t526 + t508 * t518);
t554 = t533 * t549;
t556 = t34 * t549;
t562 = t507 * (-mrSges(1,1) * t531 + t549 * t38 - t513 * t534 - t537 * t536) + (mrSges(1,1) * t549 + t38 * t531 + t513 * t554 + t537 * t556) * t504;
t567 = -t504 * t508 + t507 * t505;
t571 = t512 * t567 * t537 + t519 * t513;
t576 = t512 * t518;
t578 = t512 * t522;
t580 = g(1) * t519 - g(2) * t576 + g(3) * t578;
t589 = t537 * (-t34 * t580 - t504 * t534 - t507 * t554) - (-t504 * t536 - t507 * t556 + t533 * t580) * t513;
t594 = t512 * t567 * t513 - t537 * t519;
t601 = t531 * t504 * t513 + t549 * t513 * t507 - t580 * t537;
t604 = -t62 * t12 * t10 * t8 * t1 + t89 * t1 * t71 - t101 * m(3) * t94 - t162 * t114 * t112 * t110 * t103 + t189 * t103 * t171 - t201 * m(3) * t194 - t262 * t214 * t212 * t210 * t203 + t289 * t203 * t271 - t301 * m(3) * t294 - t362 * t314 * t312 * t310 * t303 + t389 * t303 * t371 - t401 * m(3) * t394 - t462 * t414 * t412 * t410 * t403 + t489 * t403 * t471 - t501 * m(3) * t494 - t562 * t514 * t512 * t510 * t503 + t589 * t503 * t571 - t601 * m(3) * t594 - g(1) * m(4);
t607 = -t8 * t18 + t20 * t67;
t609 = t62 * t12;
t611 = t67 * t17;
t614 = t16 * t611 + t8 * t20;
t618 = -t10 * t11 * t16 + t36 * t614;
t623 = t11 * t614 + t36 * t76;
t628 = -t110 * t120 + t122 * t167;
t630 = t162 * t114;
t632 = t167 * t119;
t635 = t110 * t122 + t118 * t632;
t639 = -t112 * t113 * t118 + t137 * t635;
t644 = t113 * t635 + t137 * t176;
t649 = -t210 * t220 + t222 * t267;
t651 = t262 * t214;
t653 = t267 * t219;
t656 = t210 * t222 + t218 * t653;
t660 = -t212 * t213 * t218 + t237 * t656;
t665 = t213 * t656 + t237 * t276;
t670 = -t310 * t320 + t322 * t367;
t672 = t362 * t314;
t674 = t367 * t319;
t677 = t310 * t322 + t318 * t674;
t681 = -t312 * t313 * t318 + t337 * t677;
t686 = t313 * t677 + t337 * t376;
t691 = -t410 * t420 + t422 * t467;
t693 = t462 * t414;
t695 = t467 * t419;
t698 = t410 * t422 + t418 * t695;
t702 = -t412 * t413 * t418 + t437 * t698;
t707 = t413 * t698 + t437 * t476;
t712 = -t510 * t520 + t522 * t567;
t714 = t562 * t514;
t716 = t567 * t519;
t719 = t510 * t522 + t518 * t716;
t723 = -t512 * t513 * t518 + t537 * t719;
t728 = t513 * t719 + t537 * t576;
t732 = t609 * t1 * t607 + t89 * t1 * t618 - t101 * m(3) * t623 + t630 * t103 * t628 + t189 * t103 * t639 - t201 * m(3) * t644 + t651 * t203 * t649 + t289 * t203 * t660 - t301 * m(3) * t665 + t672 * t303 * t670 + t389 * t303 * t681 - t401 * m(3) * t686 + t693 * t403 * t691 + t489 * t403 * t702 - t501 * m(3) * t707 + t714 * t503 * t712 + t589 * t503 * t723 - t601 * m(3) * t728 - g(2) * m(4);
t736 = t20 * t8 * t17 + t67 * t16;
t741 = t8 * t16 - t20 * t611;
t745 = t10 * t11 * t20 + t36 * t741;
t750 = t11 * t741 - t36 * t78;
t756 = t122 * t110 * t119 + t167 * t118;
t761 = t110 * t118 - t122 * t632;
t765 = t112 * t113 * t122 + t137 * t761;
t770 = t113 * t761 - t137 * t178;
t776 = t222 * t210 * t219 + t267 * t218;
t781 = t210 * t218 - t222 * t653;
t785 = t212 * t213 * t222 + t237 * t781;
t790 = t213 * t781 - t237 * t278;
t796 = t322 * t310 * t319 + t367 * t318;
t801 = t310 * t318 - t322 * t674;
t805 = t312 * t313 * t322 + t337 * t801;
t810 = t313 * t801 - t337 * t378;
t816 = t422 * t410 * t419 + t467 * t418;
t821 = t410 * t418 - t422 * t695;
t825 = t412 * t413 * t422 + t437 * t821;
t830 = t413 * t821 - t437 * t478;
t836 = t522 * t510 * t519 + t567 * t518;
t841 = t510 * t518 - t522 * t716;
t845 = t512 * t513 * t522 + t537 * t841;
t850 = t513 * t841 - t537 * t578;
t854 = t609 * t1 * t736 + t89 * t1 * t745 - t101 * m(3) * t750 + t630 * t103 * t756 + t189 * t103 * t765 - t201 * m(3) * t770 + t651 * t203 * t776 + t289 * t203 * t785 - t301 * m(3) * t790 + t672 * t303 * t796 + t389 * t303 * t805 - t401 * m(3) * t810 + t693 * t403 * t816 + t489 * t403 * t825 - t501 * m(3) * t830 + t714 * t503 * t836 + t589 * t503 * t845 - t601 * m(3) * t850 - g(3) * m(4);
t855 = sin(xP(6));
t857 = cos(xP(6));
t859 = -koppelP(1,2) * t855 + t857 * koppelP(1,1);
t860 = sin(xP(5));
t862 = cos(xP(5));
t863 = koppelP(1,3) * t862;
t864 = t860 * t859 - t863;
t865 = cos(xP(4));
t867 = sin(xP(4));
t870 = koppelP(1,1) * t855 + t857 * koppelP(1,2);
t871 = t870 * t867;
t872 = t865 * t864 - t871;
t874 = t12 * t1;
t877 = t870 * t865;
t878 = t867 * t864 + t877;
t896 = -koppelP(2,2) * t855 + t857 * koppelP(2,1);
t898 = koppelP(2,3) * t862;
t899 = t860 * t896 - t898;
t903 = koppelP(2,1) * t855 + t857 * koppelP(2,2);
t904 = t903 * t867;
t905 = t865 * t899 - t904;
t907 = t114 * t103;
t910 = t903 * t865;
t911 = t867 * t899 + t910;
t929 = -koppelP(3,2) * t855 + t857 * koppelP(3,1);
t931 = koppelP(3,3) * t862;
t932 = t860 * t929 - t931;
t936 = koppelP(3,1) * t855 + t857 * koppelP(3,2);
t937 = t936 * t867;
t938 = t865 * t932 - t937;
t940 = t214 * t203;
t943 = t936 * t865;
t944 = t867 * t932 + t943;
t962 = -koppelP(4,2) * t855 + t857 * koppelP(4,1);
t964 = koppelP(4,3) * t862;
t965 = t860 * t962 - t964;
t969 = koppelP(4,1) * t855 + t857 * koppelP(4,2);
t970 = t969 * t867;
t971 = t865 * t965 - t970;
t973 = t314 * t303;
t976 = t969 * t865;
t977 = t867 * t965 + t976;
t995 = -koppelP(5,2) * t855 + t857 * koppelP(5,1);
t997 = koppelP(5,3) * t862;
t998 = t860 * t995 - t997;
t1002 = koppelP(5,1) * t855 + t857 * koppelP(5,2);
t1003 = t1002 * t867;
t1004 = t865 * t998 - t1003;
t1006 = t414 * t403;
t1009 = t1002 * t865;
t1010 = t867 * t998 + t1009;
t1028 = -koppelP(6,2) * t855 + t857 * koppelP(6,1);
t1030 = koppelP(6,3) * t862;
t1031 = t860 * t1028 - t1030;
t1035 = koppelP(6,1) * t855 + t857 * koppelP(6,2);
t1036 = t1035 * t867;
t1037 = t865 * t1031 - t1036;
t1039 = t514 * t503;
t1042 = t1035 * t865;
t1043 = t867 * t1031 + t1042;
t1059 = t860 * g(2);
t1061 = g(3) * mrSges(4,2);
t1064 = mrSges(4,2) * t855;
t1072 = g(3) * t860;
t1079 = t855 * g(2);
t1085 = t62 * (t874 * t607 * t872 + t874 * t736 * t878) + t89 * (t1 * t618 * t872 + t1 * t745 * t878) - t101 * m(3) * (t623 * t872 + t750 * t878) + t162 * (t907 * t628 * t905 + t907 * t756 * t911) + t189 * (t103 * t639 * t905 + t103 * t765 * t911) - t201 * m(3) * (t644 * t905 + t770 * t911) + t262 * (t940 * t649 * t938 + t940 * t776 * t944) + t289 * (t203 * t660 * t938 + t203 * t785 * t944) - t301 * m(3) * (t665 * t938 + t790 * t944) + t362 * (t973 * t670 * t971 + t973 * t796 * t977) + t389 * (t303 * t681 * t971 + t303 * t805 * t977) - t401 * m(3) * (t686 * t971 + t810 * t977) + t462 * (t1006 * t691 * t1004 + t1006 * t816 * t1010) + t489 * (t403 * t702 * t1004 + t403 * t825 * t1010) - t501 * m(3) * (t707 * t1004 + t830 * t1010) + t562 * (t1039 * t712 * t1037 + t1039 * t836 * t1043) + t589 * (t503 * t723 * t1037 + t503 * t845 * t1043) - t601 * m(3) * (t728 * t1037 + t850 * t1043) - t865 * (t857 * (mrSges(4,1) * t1059 + t1061) - t1064 * t1059 - mrSges(4,3) * t862 * g(2) + t855 * mrSges(4,1) * g(3)) + (t857 * (g(2) * mrSges(4,2) - mrSges(4,1) * t1072) + t1061 * t855 * t860 + mrSges(4,1) * t1079 + mrSges(4,3) * g(3) * t862) * t867;
t1087 = -t860 * t859 + t863;
t1089 = t865 * t1087 + t871;
t1092 = t12 * t10 * t8;
t1094 = t857 * t862;
t1096 = t855 * t862;
t1099 = -koppelP(1,1) * t1094 + koppelP(1,2) * t1096 - koppelP(1,3) * t860;
t1116 = -t860 * t896 + t898;
t1118 = t865 * t1116 + t904;
t1121 = t114 * t112 * t110;
t1126 = -koppelP(2,1) * t1094 + koppelP(2,2) * t1096 - koppelP(2,3) * t860;
t1143 = -t860 * t929 + t931;
t1145 = t865 * t1143 + t937;
t1148 = t214 * t212 * t210;
t1153 = -koppelP(3,1) * t1094 + koppelP(3,2) * t1096 - koppelP(3,3) * t860;
t1170 = -t860 * t962 + t964;
t1172 = t865 * t1170 + t970;
t1175 = t314 * t312 * t310;
t1180 = -koppelP(4,1) * t1094 + koppelP(4,2) * t1096 - koppelP(4,3) * t860;
t1191 = t62 * (-t1092 * t1 * t1089 + t874 * t736 * t1099) + t89 * (t1 * t71 * t1089 + t1 * t745 * t1099) - t101 * m(3) * (t94 * t1089 + t750 * t1099) + t162 * (-t1121 * t103 * t1118 + t907 * t756 * t1126) + t189 * (t103 * t171 * t1118 + t103 * t765 * t1126) - t201 * m(3) * (t194 * t1118 + t770 * t1126) + t262 * (-t1148 * t203 * t1145 + t940 * t776 * t1153) + t289 * (t203 * t271 * t1145 + t203 * t785 * t1153) - t301 * m(3) * (t294 * t1145 + t790 * t1153) + t362 * (-t1175 * t303 * t1172 + t973 * t796 * t1180) + t389 * (t303 * t371 * t1172 + t303 * t805 * t1180);
t1198 = -t860 * t995 + t997;
t1200 = t865 * t1198 + t1003;
t1203 = t414 * t412 * t410;
t1208 = -koppelP(5,1) * t1094 + koppelP(5,2) * t1096 - koppelP(5,3) * t860;
t1225 = -t860 * t1028 + t1030;
t1227 = t865 * t1225 + t1036;
t1230 = t514 * t512 * t510;
t1235 = -koppelP(6,1) * t1094 + koppelP(6,2) * t1096 - koppelP(6,3) * t860;
t1253 = mrSges(4,1) * t857 - t1064;
t1264 = g(1) * mrSges(4,1);
t1267 = -t401 * m(3) * (t394 * t1172 + t810 * t1180) + t462 * (t1006 * t816 * t1208 - t1203 * t403 * t1200) + t489 * (t403 * t471 * t1200 + t403 * t825 * t1208) - t501 * m(3) * (t494 * t1200 + t830 * t1208) + t562 * (t1039 * t836 * t1235 - t1230 * t503 * t1227) + t589 * (t503 * t571 * t1227 + t503 * t845 * t1235) - t601 * m(3) * (t594 * t1227 + t850 * t1235) - t865 * (-t1253 * t860 + mrSges(4,3) * t862) * g(1) + t862 * t1253 * g(3) - g(1) * mrSges(4,2) * t867 * t857 - t1264 * t855 * t867 + mrSges(4,3) * t1072;
t1270 = t867 * t1087 - t877;
t1289 = t867 * t1116 - t910;
t1308 = t867 * t1143 - t943;
t1327 = t867 * t1170 - t976;
t1340 = t62 * (-t1092 * t1 * t1270 - t874 * t607 * t1099) + t89 * (-t1 * t618 * t1099 + t1 * t71 * t1270) - t101 * m(3) * (-t623 * t1099 + t94 * t1270) + t162 * (-t1121 * t103 * t1289 - t907 * t628 * t1126) + t189 * (-t103 * t639 * t1126 + t103 * t171 * t1289) - t201 * m(3) * (-t644 * t1126 + t194 * t1289) + t262 * (-t1148 * t203 * t1308 - t940 * t649 * t1153) + t289 * (-t203 * t660 * t1153 + t203 * t271 * t1308) - t301 * m(3) * (-t665 * t1153 + t294 * t1308) + t362 * (-t1175 * t303 * t1327 - t973 * t670 * t1180) + t389 * (-t303 * t681 * t1180 + t303 * t371 * t1327);
t1347 = t867 * t1198 - t1009;
t1366 = t867 * t1225 - t1042;
t1391 = t860 * t867;
t1403 = -t401 * m(3) * (-t686 * t1180 + t394 * t1327) + t462 * (-t1006 * t691 * t1208 - t1203 * t403 * t1347) + t489 * (-t403 * t702 * t1208 + t403 * t471 * t1347) - t501 * m(3) * (-t707 * t1208 + t494 * t1347) + t562 * (-t1039 * t712 * t1235 - t1230 * t503 * t1366) + t589 * (-t503 * t723 * t1235 + t503 * t571 * t1366) - t601 * m(3) * (-t728 * t1235 + t594 * t1366) - t862 * (mrSges(4,3) * t867 * g(1) + t857 * mrSges(4,1) * g(2) - mrSges(4,2) * t1079) + t857 * (mrSges(4,1) * t1391 + t865 * mrSges(4,2)) * g(1) - mrSges(4,2) * g(1) * t855 * t1391 + t1264 * t855 * t865 - mrSges(4,3) * t1059;
unknown(1,1) = t604;
unknown(2,1) = t732;
unknown(3,1) = t854;
unknown(4,1) = t1085;
unknown(5,1) = t1191 + t1267;
unknown(6,1) = t1340 + t1403;
taugX  = unknown;
