% Calculate inertia matrix for parallel robot
% P3PRRRR8V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V1G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:49:56
% EndTime: 2020-08-06 16:49:57
% DurationCPUTime: 0.99s
% Computational Cost: add. (2133->169), mult. (5223->350), div. (360->7), fcn. (5856->22), ass. (0->178)
t460 = cos(qJ(3,1));
t440 = 0.1e1 / t460;
t455 = sin(qJ(2,1));
t461 = cos(qJ(2,1));
t533 = pkin(2) * t460;
t418 = -t461 * pkin(5) + t455 * t533;
t442 = sin(pkin(3));
t444 = cos(pkin(3));
t454 = sin(qJ(3,1));
t536 = pkin(2) * t454;
t466 = 0.1e1 / (t418 * t442 + t444 * t536);
t523 = t466 * t440;
t458 = cos(qJ(3,2));
t439 = 0.1e1 / t458;
t453 = sin(qJ(2,2));
t459 = cos(qJ(2,2));
t534 = pkin(2) * t458;
t417 = -t459 * pkin(5) + t453 * t534;
t452 = sin(qJ(3,2));
t537 = pkin(2) * t452;
t467 = 0.1e1 / (t417 * t442 + t444 * t537);
t524 = t467 * t439;
t456 = cos(qJ(3,3));
t438 = 0.1e1 / t456;
t451 = sin(qJ(2,3));
t457 = cos(qJ(2,3));
t535 = pkin(2) * t456;
t416 = -t457 * pkin(5) + t451 * t535;
t450 = sin(qJ(3,3));
t538 = pkin(2) * t450;
t468 = 0.1e1 / (t416 * t442 + t444 * t538);
t525 = t468 * t438;
t539 = 2 * Ifges(3,4);
t532 = Ifges(3,1) + Ifges(2,3);
t445 = legFrame(3,3);
t428 = sin(t445);
t431 = cos(t445);
t441 = sin(pkin(6));
t443 = cos(pkin(6));
t401 = -t441 * t428 + t431 * t443;
t507 = t444 * t451;
t407 = t441 * t457 + t443 * t507;
t410 = -t441 * t507 + t443 * t457;
t510 = t442 * t456;
t377 = (-t407 * t431 - t428 * t410) * t450 - t401 * t510;
t531 = t377 * t438;
t446 = legFrame(2,3);
t429 = sin(t446);
t432 = cos(t446);
t402 = -t441 * t429 + t432 * t443;
t506 = t444 * t453;
t408 = t441 * t459 + t443 * t506;
t411 = -t441 * t506 + t443 * t459;
t509 = t442 * t458;
t378 = (-t408 * t432 - t429 * t411) * t452 - t402 * t509;
t530 = t378 * t439;
t447 = legFrame(1,3);
t430 = sin(t447);
t433 = cos(t447);
t403 = -t441 * t430 + t433 * t443;
t505 = t444 * t455;
t409 = t441 * t461 + t443 * t505;
t412 = -t441 * t505 + t443 * t461;
t508 = t442 * t460;
t379 = (-t409 * t433 - t430 * t412) * t454 - t403 * t508;
t529 = t379 * t440;
t404 = t443 * t428 + t431 * t441;
t380 = (-t428 * t407 + t410 * t431) * t450 - t404 * t510;
t528 = t380 * t438;
t405 = t443 * t429 + t432 * t441;
t381 = (-t429 * t408 + t411 * t432) * t452 - t405 * t509;
t527 = t381 * t439;
t406 = t443 * t430 + t433 * t441;
t382 = (-t430 * t409 + t412 * t433) * t454 - t406 * t508;
t526 = t382 * t440;
t437 = m(1) + m(2) + m(3);
t522 = t468 * t437;
t521 = t467 * t437;
t520 = t466 * t437;
t449 = mrSges(2,2) - mrSges(3,3);
t483 = t456 * mrSges(3,1) - t450 * mrSges(3,2);
t519 = ((mrSges(2,1) + t483) * t457 - t451 * t449) * t442;
t482 = t458 * mrSges(3,1) - t452 * mrSges(3,2);
t518 = ((mrSges(2,1) + t482) * t459 - t453 * t449) * t442;
t481 = t460 * mrSges(3,1) - t454 * mrSges(3,2);
t517 = ((mrSges(2,1) + t481) * t461 - t455 * t449) * t442;
t448 = -Ifges(3,1) + Ifges(3,2);
t516 = ((t448 * t456 + t450 * t539) * t456 + t532) * t438;
t515 = ((t448 * t458 + t452 * t539) * t458 + t532) * t439;
t514 = ((t448 * t460 + t454 * t539) * t460 + t532) * t440;
t422 = Ifges(3,5) * t450 + Ifges(3,6) * t456;
t513 = t422 * t438;
t423 = Ifges(3,5) * t452 + Ifges(3,6) * t458;
t512 = t423 * t439;
t424 = Ifges(3,5) * t454 + Ifges(3,6) * t460;
t511 = t424 * t440;
t504 = t451 * t401;
t503 = t451 * t404;
t502 = t453 * t402;
t501 = t453 * t405;
t500 = t455 * t403;
t499 = t455 * t406;
t498 = t457 * t401;
t497 = t457 * t404;
t496 = t459 * t402;
t495 = t459 * t405;
t494 = t461 * t403;
t493 = t461 * t406;
t371 = -(t444 * t498 - t503) * t535 - pkin(5) * (t444 * t504 + t497);
t492 = t371 * t525;
t372 = -(t444 * t497 + t504) * t535 - (t444 * t503 - t498) * pkin(5);
t491 = t372 * t525;
t373 = -(t444 * t496 - t501) * t534 - pkin(5) * (t444 * t502 + t495);
t490 = t373 * t524;
t374 = -(t444 * t495 + t502) * t534 - (t444 * t501 - t496) * pkin(5);
t489 = t374 * t524;
t375 = -(t444 * t494 - t499) * t533 - pkin(5) * (t444 * t500 + t493);
t488 = t375 * t523;
t376 = -(t444 * t493 + t500) * t533 - (t444 * t499 - t494) * pkin(5);
t487 = t376 * t523;
t462 = 0.1e1 / pkin(2);
t486 = t462 * t525;
t485 = t462 * t524;
t484 = t462 * t523;
t419 = pkin(5) * t451 + t457 * t535;
t465 = -t416 * t444 + t442 * t538;
t383 = t419 * t443 + t465 * t441;
t386 = t441 * t419 - t465 * t443;
t365 = t383 * t431 - t428 * t386;
t471 = t519 * t525;
t389 = t483 * t444 - (t450 * mrSges(3,1) + t456 * mrSges(3,2)) * t442 * t451;
t477 = t389 * t486;
t329 = t365 * t522 + t371 * t477 + t377 * t471;
t420 = pkin(5) * t453 + t459 * t534;
t464 = -t417 * t444 + t442 * t537;
t384 = t420 * t443 + t464 * t441;
t387 = t441 * t420 - t464 * t443;
t366 = t384 * t432 - t429 * t387;
t470 = t518 * t524;
t390 = t482 * t444 - (t452 * mrSges(3,1) + t458 * mrSges(3,2)) * t442 * t453;
t476 = t390 * t485;
t330 = t366 * t521 + t373 * t476 + t378 * t470;
t421 = pkin(5) * t455 + t461 * t533;
t463 = -t418 * t444 + t442 * t536;
t385 = t421 * t443 + t463 * t441;
t388 = t441 * t421 - t463 * t443;
t367 = t385 * t433 - t430 * t388;
t469 = t517 * t523;
t391 = t481 * t444 - (t454 * mrSges(3,1) + t460 * mrSges(3,2)) * t442 * t455;
t475 = t391 * t484;
t331 = t367 * t520 + t375 * t475 + t379 * t469;
t368 = t383 * t428 + t386 * t431;
t332 = t368 * t522 + t372 * t477 + t380 * t471;
t369 = t384 * t429 + t387 * t432;
t333 = t369 * t521 + t374 * t476 + t381 * t470;
t370 = t385 * t430 + t388 * t433;
t334 = t370 * t520 + t376 * t475 + t382 * t469;
t480 = Ifges(3,3) * t486;
t479 = Ifges(3,3) * t485;
t478 = Ifges(3,3) * t484;
t474 = t422 * t486;
t473 = t423 * t485;
t472 = t424 * t484;
t346 = t376 * t478 + (t370 * t391 + t382 * t511) * t466;
t345 = t374 * t479 + (t369 * t390 + t381 * t512) * t467;
t344 = t372 * t480 + (t368 * t389 + t380 * t513) * t468;
t343 = t375 * t478 + (t367 * t391 + t379 * t511) * t466;
t342 = t373 * t479 + (t366 * t390 + t378 * t512) * t467;
t341 = t371 * t480 + (t365 * t389 + t377 * t513) * t468;
t340 = t376 * t472 + (t370 * t517 + t382 * t514) * t466;
t339 = t374 * t473 + (t369 * t518 + t381 * t515) * t467;
t338 = t372 * t474 + (t368 * t519 + t380 * t516) * t468;
t337 = t375 * t472 + (t367 * t517 + t379 * t514) * t466;
t336 = t373 * t473 + (t366 * t518 + t378 * t515) * t467;
t335 = t371 * t474 + (t365 * t519 + t377 * t516) * t468;
t328 = t334 + t333 + t332;
t327 = t331 + t330 + t329;
t1 = [m(4) + (t331 * t367 + t337 * t529) * t466 + (t330 * t366 + t336 * t530) * t467 + (t329 * t365 + t335 * t531) * t468 + (t341 * t492 + t342 * t490 + t343 * t488) * t462, (t331 * t370 + t337 * t526) * t466 + (t330 * t369 + t336 * t527) * t467 + (t329 * t368 + t335 * t528) * t468 + (t341 * t491 + t342 * t489 + t343 * t487) * t462, t327; (t334 * t367 + t340 * t529) * t466 + (t333 * t366 + t339 * t530) * t467 + (t332 * t365 + t338 * t531) * t468 + (t344 * t492 + t345 * t490 + t346 * t488) * t462, m(4) + (t334 * t370 + t340 * t526) * t466 + (t333 * t369 + t339 * t527) * t467 + (t332 * t368 + t338 * t528) * t468 + (t344 * t491 + t345 * t489 + t346 * t487) * t462, t328; t327, t328, (3 * m(1)) + (3 * m(2)) + (3 * m(3)) + m(4);];
MX  = t1;
