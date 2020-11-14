% Calculate inertia matrix for parallel robot
% P3RRPRR12V1G1A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
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
% Datum: 2020-08-06 19:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V1G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:01:46
% EndTime: 2020-08-06 19:01:47
% DurationCPUTime: 0.76s
% Computational Cost: add. (2013->203), mult. (2907->357), div. (330->6), fcn. (2355->18), ass. (0->177)
t465 = (qJ(3,3) ^ 2);
t534 = m(3) * t465;
t467 = (qJ(3,2) ^ 2);
t533 = m(3) * t467;
t469 = (qJ(3,1) ^ 2);
t532 = m(3) * t469;
t451 = sin(qJ(1,3));
t531 = t451 * pkin(4);
t453 = sin(qJ(1,2));
t530 = t453 * pkin(4);
t455 = sin(qJ(1,1));
t529 = t455 * pkin(4);
t528 = Ifges(2,1) + Ifges(3,1);
t527 = Ifges(2,6) - Ifges(3,6);
t450 = sin(qJ(2,3));
t526 = mrSges(3,2) * t450;
t452 = sin(qJ(2,2));
t525 = mrSges(3,2) * t452;
t454 = sin(qJ(2,1));
t524 = mrSges(3,2) * t454;
t523 = (mrSges(3,3) * qJ(3,1));
t522 = (mrSges(3,3) * qJ(3,2));
t521 = (mrSges(3,3) * qJ(3,3));
t456 = cos(qJ(2,3));
t463 = pkin(1) + pkin(2);
t499 = t463 * t456;
t508 = t450 * qJ(3,3);
t429 = t499 + t508;
t457 = cos(qJ(1,3));
t440 = t457 * pkin(4);
t399 = t429 * t451 + t440;
t402 = t429 * t457 - t531;
t447 = legFrame(3,3);
t433 = sin(t447);
t436 = cos(t447);
t381 = -t433 * t399 + t436 * t402;
t466 = 0.1e1 / qJ(3,3);
t520 = t381 * t466;
t458 = cos(qJ(2,2));
t498 = t463 * t458;
t506 = t452 * qJ(3,2);
t430 = t498 + t506;
t459 = cos(qJ(1,2));
t441 = t459 * pkin(4);
t400 = t430 * t453 + t441;
t403 = t430 * t459 - t530;
t448 = legFrame(2,3);
t434 = sin(t448);
t437 = cos(t448);
t382 = -t434 * t400 + t437 * t403;
t468 = 0.1e1 / qJ(3,2);
t519 = t382 * t468;
t460 = cos(qJ(2,1));
t497 = t463 * t460;
t504 = t454 * qJ(3,1);
t431 = t497 + t504;
t461 = cos(qJ(1,1));
t442 = t461 * pkin(4);
t401 = t431 * t455 + t442;
t404 = t431 * t461 - t529;
t449 = legFrame(1,3);
t435 = sin(t449);
t438 = cos(t449);
t383 = -t435 * t401 + t438 * t404;
t470 = 0.1e1 / qJ(3,1);
t518 = t383 * t470;
t384 = t399 * t436 + t433 * t402;
t517 = t384 * t466;
t385 = t400 * t437 + t434 * t403;
t516 = t385 * t468;
t386 = t401 * t438 + t435 * t404;
t515 = t386 * t470;
t426 = -t456 * qJ(3,3) + t463 * t450;
t439 = m(3) * pkin(1) + mrSges(3,1);
t396 = (m(3) * t426 - t450 * t439) * t466;
t514 = t396 * t466;
t427 = -t458 * qJ(3,2) + t463 * t452;
t397 = (m(3) * t427 - t452 * t439) * t468;
t513 = t397 * t468;
t428 = -t460 * qJ(3,1) + t463 * t454;
t398 = (m(3) * t428 - t454 * t439) * t470;
t512 = t398 * t470;
t511 = t439 * t466;
t510 = t439 * t468;
t509 = t439 * t470;
t507 = t450 * t466;
t505 = t452 * t468;
t503 = t454 * t470;
t502 = t456 * t466;
t501 = t458 * t468;
t500 = t460 * t470;
t496 = Ifges(1,3) + t528;
t495 = pkin(1) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t446 = 2 * mrSges(3,1) * pkin(1);
t494 = Ifges(3,2) + Ifges(2,3) + t446;
t493 = mrSges(3,2) * t507;
t492 = mrSges(3,2) * t505;
t491 = mrSges(3,2) * t503;
t408 = -t433 * t451 + t436 * t457;
t417 = t451 * t508 + t440;
t418 = t457 * t508 - t531;
t375 = t408 * t499 - t433 * t417 + t418 * t436;
t490 = t375 * t502;
t409 = t433 * t457 + t436 * t451;
t376 = t409 * t499 + t417 * t436 + t433 * t418;
t489 = t376 * t502;
t410 = -t434 * t453 + t437 * t459;
t419 = t453 * t506 + t441;
t420 = t459 * t506 - t530;
t377 = t410 * t498 - t434 * t419 + t420 * t437;
t488 = t377 * t501;
t411 = t434 * t459 + t437 * t453;
t378 = t411 * t498 + t419 * t437 + t434 * t420;
t487 = t378 * t501;
t412 = -t435 * t455 + t438 * t461;
t421 = t455 * t504 + t442;
t422 = t461 * t504 - t529;
t379 = t412 * t497 - t435 * t421 + t422 * t438;
t486 = t379 * t500;
t413 = t435 * t461 + t438 * t455;
t380 = t413 * t497 + t421 * t438 + t435 * t422;
t485 = t380 * t500;
t443 = 2 * t521;
t471 = pkin(1) ^ 2;
t414 = (t465 + t471) * m(3) + t443 + t494;
t393 = (t414 * t450 - t426 * t439) * t466;
t484 = t393 * t502;
t444 = 2 * t522;
t415 = (t467 + t471) * m(3) + t444 + t494;
t394 = (t415 * t452 - t427 * t439) * t468;
t483 = t394 * t501;
t445 = 2 * t523;
t416 = (t469 + t471) * m(3) + t445 + t494;
t395 = (t416 * t454 - t428 * t439) * t470;
t482 = t395 * t500;
t432 = pkin(1) * mrSges(3,2) - Ifges(3,4) - Ifges(2,5);
t405 = (mrSges(3,2) * qJ(3,3) + t527) * t456 - t450 * t432;
t481 = t405 * t502;
t406 = (mrSges(3,2) * qJ(3,2) + t527) * t458 - t452 * t432;
t480 = t406 * t501;
t407 = (qJ(3,1) * mrSges(3,2) + t527) * t460 - t454 * t432;
t479 = t407 * t500;
t478 = t414 * t502;
t477 = t415 * t501;
t476 = t416 * t500;
t475 = t439 * t502;
t474 = t439 * t501;
t473 = t439 * t500;
t472 = m(3) * t471 + Ifges(2,2) + Ifges(3,3) + t446 - t528;
t425 = 0.1e1 / t431;
t424 = 0.1e1 / t430;
t423 = 0.1e1 / t429;
t392 = (mrSges(3,2) * t428 + t407) * t503;
t391 = (mrSges(3,2) * t427 + t406) * t505;
t390 = (mrSges(3,2) * t426 + t405) * t507;
t389 = t532 + t445 + (0.2e1 * t454 * (t439 * qJ(3,1) + t495) + (t472 - (2 * t523) - t532) * t460) * t460 + t496;
t388 = t533 + t444 + (0.2e1 * t452 * (t439 * qJ(3,2) + t495) + (t472 - (2 * t522) - t533) * t458) * t458 + t496;
t387 = t534 + t443 + (0.2e1 * t450 * (t439 * qJ(3,3) + t495) + (t472 - (2 * t521) - t534) * t456) * t456 + t496;
t374 = m(3) * t518 + (-t379 * t473 - t413 * t524) * t425;
t373 = m(3) * t515 + (-t380 * t473 + t412 * t524) * t425;
t372 = m(3) * t519 + (-t377 * t474 - t411 * t525) * t424;
t371 = m(3) * t516 + (-t378 * t474 + t410 * t525) * t424;
t370 = m(3) * t520 + (-t375 * t475 - t409 * t526) * t423;
t369 = m(3) * t517 + (-t376 * t475 + t408 * t526) * t423;
t368 = -t383 * t509 + (t379 * t476 - t407 * t413) * t425;
t367 = -t386 * t509 + (t380 * t476 + t407 * t412) * t425;
t366 = -t382 * t510 + (t377 * t477 - t406 * t411) * t424;
t365 = -t385 * t510 + (t378 * t477 + t406 * t410) * t424;
t364 = -t381 * t511 + (t375 * t478 - t405 * t409) * t423;
t363 = -t384 * t511 + (t376 * t478 + t405 * t408) * t423;
t362 = t383 * t491 + (t379 * t479 - t389 * t413) * t425;
t361 = t386 * t491 + (t380 * t479 + t389 * t412) * t425;
t360 = t382 * t492 + (t377 * t480 - t388 * t411) * t424;
t359 = t385 * t492 + (t378 * t480 + t388 * t410) * t424;
t358 = t381 * t493 + (t375 * t481 - t387 * t409) * t423;
t357 = t384 * t493 + (t376 * t481 + t387 * t408) * t423;
t1 = [t370 * t520 + t372 * t519 + t374 * t518 + m(4) + (-t362 * t413 + t368 * t486) * t425 + (-t360 * t411 + t366 * t488) * t424 + (-t358 * t409 + t364 * t490) * t423, t370 * t517 + t372 * t516 + t374 * t515 + (t362 * t412 + t368 * t485) * t425 + (t360 * t410 + t366 * t487) * t424 + (t358 * t408 + t364 * t489) * t423, (t368 * t454 + t374 * t428) * t470 + (t366 * t452 + t372 * t427) * t468 + (t364 * t450 + t370 * t426) * t466; t369 * t520 + t371 * t519 + t373 * t518 + (-t361 * t413 + t367 * t486) * t425 + (-t359 * t411 + t365 * t488) * t424 + (-t357 * t409 + t363 * t490) * t423, t369 * t517 + t371 * t516 + t373 * t515 + m(4) + (t361 * t412 + t367 * t485) * t425 + (t359 * t410 + t365 * t487) * t424 + (t357 * t408 + t363 * t489) * t423, (t367 * t454 + t373 * t428) * t470 + (t365 * t452 + t371 * t427) * t468 + (t363 * t450 + t369 * t426) * t466; t381 * t514 + t382 * t513 + t383 * t512 + (t379 * t482 - t392 * t413) * t425 + (t377 * t483 - t391 * t411) * t424 + (t375 * t484 - t390 * t409) * t423, t384 * t514 + t385 * t513 + t386 * t512 + (t380 * t482 + t392 * t412) * t425 + (t378 * t483 + t391 * t410) * t424 + (t376 * t484 + t390 * t408) * t423, m(4) + (t395 * t454 + t398 * t428) * t470 + (t394 * t452 + t397 * t427) * t468 + (t393 * t450 + t396 * t426) * t466;];
MX  = t1;
