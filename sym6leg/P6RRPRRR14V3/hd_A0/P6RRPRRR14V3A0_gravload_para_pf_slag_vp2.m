% Calculate Gravitation load for parallel robot
% P6RRPRRR14V3G1P1A0
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
% Datum: 2019-04-15 09:53
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P6RRPRRR14V3G1P1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(1,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G1P1A0_gravload_para_pf_slag_vp2: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G1P1A0_gravload_para_pf_slag_vp2: qJ has to be [3x6] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G1P1A0_gravload_para_pf_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6RRPRRR14V3G1P1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6RRPRRR14V3G1P1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P6RRPRRR14V3G1P1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G1P1A0_gravload_para_pf_slag_vp2: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G1P1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-15 09:51:27
% EndTime: 2019-04-15 09:51:31
% DurationCPUTime: 4.62s
% Computational Cost: add. (2141->281), mult. (3378->512), div. (126->12), fcn. (2783->42), ass. (0->233)
t3493 = xP(5);
t3449 = sin(t3493);
t3452 = cos(t3493);
t3504 = koppelP(6,3);
t3492 = xP(6);
t3448 = sin(t3492);
t3451 = cos(t3492);
t3510 = koppelP(6,2);
t3516 = koppelP(6,1);
t3620 = -t3448 * t3510 + t3451 * t3516;
t3375 = -t3449 * t3620 + t3452 * t3504;
t3406 = t3448 * t3516 + t3451 * t3510;
t3494 = xP(4);
t3450 = sin(t3494);
t3453 = cos(t3494);
t3343 = t3375 * t3450 - t3406 * t3453;
t3339 = t3375 * t3453 + t3406 * t3450;
t3505 = koppelP(5,3);
t3511 = koppelP(5,2);
t3517 = koppelP(5,1);
t3621 = -t3448 * t3511 + t3451 * t3517;
t3377 = -t3449 * t3621 + t3452 * t3505;
t3407 = t3448 * t3517 + t3451 * t3511;
t3346 = t3377 * t3450 - t3407 * t3453;
t3340 = t3377 * t3453 + t3407 * t3450;
t3506 = koppelP(4,3);
t3512 = koppelP(4,2);
t3518 = koppelP(4,1);
t3622 = -t3448 * t3512 + t3451 * t3518;
t3379 = -t3449 * t3622 + t3452 * t3506;
t3408 = t3448 * t3518 + t3451 * t3512;
t3349 = t3379 * t3450 - t3408 * t3453;
t3341 = t3379 * t3453 + t3408 * t3450;
t3507 = koppelP(3,3);
t3513 = koppelP(3,2);
t3519 = koppelP(3,1);
t3623 = -t3448 * t3513 + t3451 * t3519;
t3381 = -t3449 * t3623 + t3452 * t3507;
t3409 = t3448 * t3519 + t3451 * t3513;
t3352 = t3381 * t3450 - t3409 * t3453;
t3342 = t3381 * t3453 + t3409 * t3450;
t3508 = koppelP(2,3);
t3514 = koppelP(2,2);
t3520 = koppelP(2,1);
t3624 = -t3448 * t3514 + t3451 * t3520;
t3383 = -t3449 * t3624 + t3452 * t3508;
t3410 = t3448 * t3520 + t3451 * t3514;
t3355 = t3383 * t3450 - t3410 * t3453;
t3356 = t3383 * t3453 + t3410 * t3450;
t3509 = koppelP(1,3);
t3515 = koppelP(1,2);
t3521 = koppelP(1,1);
t3625 = -t3448 * t3515 + t3451 * t3521;
t3385 = -t3449 * t3625 + t3452 * t3509;
t3411 = t3448 * t3521 + t3451 * t3515;
t3359 = t3385 * t3450 - t3411 * t3453;
t3360 = t3385 * t3453 + t3411 * t3450;
t3496 = mrSges(4,2);
t3497 = mrSges(4,1);
t3528 = t3448 * t3496 - t3451 * t3497;
t3495 = mrSges(4,3);
t3554 = t3452 * t3495;
t3619 = t3449 * t3528 + t3554;
t3626 = -t3449 * t3495 + t3452 * t3528;
t3618 = mrSges(3,3) - mrSges(2,2);
t3479 = mrSges(2,1) + mrSges(3,1);
t3617 = g(3) * t3479;
t3461 = legFrame(6,3);
t3436 = sin(t3461);
t3442 = cos(t3461);
t3412 = -g(1) * t3436 + g(2) * t3442;
t3418 = g(1) * t3442 + g(2) * t3436;
t3454 = mrSges(1,2) - mrSges(3,2) - mrSges(2,3);
t3467 = sin(qJ(2,6));
t3468 = sin(qJ(1,6));
t3474 = cos(qJ(1,6));
t3430 = m(3) * qJ(3,6) + t3618;
t3473 = cos(qJ(2,6));
t3527 = t3430 * t3467 + t3473 * t3479 + mrSges(1,1);
t3616 = ((t3454 * t3474 + t3468 * t3527) * t3418 + (t3454 * t3468 - t3474 * t3527) * t3412) / t3467;
t3462 = legFrame(5,3);
t3437 = sin(t3462);
t3443 = cos(t3462);
t3413 = -g(1) * t3437 + g(2) * t3443;
t3419 = g(1) * t3443 + g(2) * t3437;
t3469 = sin(qJ(2,5));
t3470 = sin(qJ(1,5));
t3476 = cos(qJ(1,5));
t3431 = m(3) * qJ(3,5) + t3618;
t3475 = cos(qJ(2,5));
t3526 = t3431 * t3469 + t3475 * t3479 + mrSges(1,1);
t3615 = ((t3454 * t3476 + t3470 * t3526) * t3419 + (t3454 * t3470 - t3476 * t3526) * t3413) / t3469;
t3463 = legFrame(4,3);
t3438 = sin(t3463);
t3444 = cos(t3463);
t3414 = -g(1) * t3438 + g(2) * t3444;
t3420 = g(1) * t3444 + g(2) * t3438;
t3471 = sin(qJ(2,4));
t3472 = sin(qJ(1,4));
t3478 = cos(qJ(1,4));
t3432 = m(3) * qJ(3,4) + t3618;
t3477 = cos(qJ(2,4));
t3525 = t3432 * t3471 + t3477 * t3479 + mrSges(1,1);
t3614 = ((t3454 * t3478 + t3472 * t3525) * t3420 + (t3454 * t3472 - t3478 * t3525) * t3414) / t3471;
t3464 = legFrame(3,3);
t3439 = sin(t3464);
t3445 = cos(t3464);
t3415 = -g(1) * t3439 + g(2) * t3445;
t3421 = g(1) * t3445 + g(2) * t3439;
t3480 = sin(qJ(2,3));
t3481 = sin(qJ(1,3));
t3487 = cos(qJ(1,3));
t3433 = m(3) * qJ(3,3) + t3618;
t3486 = cos(qJ(2,3));
t3524 = t3433 * t3480 + t3479 * t3486 + mrSges(1,1);
t3613 = ((t3454 * t3487 + t3481 * t3524) * t3421 + (t3454 * t3481 - t3487 * t3524) * t3415) / t3480;
t3465 = legFrame(2,3);
t3440 = sin(t3465);
t3446 = cos(t3465);
t3416 = -g(1) * t3440 + g(2) * t3446;
t3422 = g(1) * t3446 + g(2) * t3440;
t3482 = sin(qJ(2,2));
t3483 = sin(qJ(1,2));
t3489 = cos(qJ(1,2));
t3434 = m(3) * qJ(3,2) + t3618;
t3488 = cos(qJ(2,2));
t3523 = t3434 * t3482 + t3479 * t3488 + mrSges(1,1);
t3612 = ((t3454 * t3489 + t3483 * t3523) * t3422 + (t3454 * t3483 - t3489 * t3523) * t3416) / t3482;
t3466 = legFrame(1,3);
t3441 = sin(t3466);
t3447 = cos(t3466);
t3417 = -g(1) * t3441 + g(2) * t3447;
t3423 = g(1) * t3447 + g(2) * t3441;
t3484 = sin(qJ(2,1));
t3485 = sin(qJ(1,1));
t3491 = cos(qJ(1,1));
t3435 = m(3) * qJ(3,1) + t3618;
t3490 = cos(qJ(2,1));
t3522 = t3435 * t3484 + t3479 * t3490 + mrSges(1,1);
t3611 = ((t3454 * t3491 + t3485 * t3522) * t3423 + (t3454 * t3485 - t3491 * t3522) * t3417) / t3484;
t3535 = t3412 * t3468 + t3418 * t3474;
t3333 = (-t3430 * t3535 - t3617) * t3473 - (g(3) * t3430 - t3479 * t3535) * t3467;
t3498 = 0.1e1 / qJ(3,6);
t3610 = t3333 * t3498;
t3534 = t3413 * t3470 + t3419 * t3476;
t3334 = (-t3431 * t3534 - t3617) * t3475 - (g(3) * t3431 - t3479 * t3534) * t3469;
t3499 = 0.1e1 / qJ(3,5);
t3609 = t3334 * t3499;
t3533 = t3414 * t3472 + t3420 * t3478;
t3335 = (-t3432 * t3533 - t3617) * t3477 - (g(3) * t3432 - t3479 * t3533) * t3471;
t3500 = 0.1e1 / qJ(3,4);
t3608 = t3335 * t3500;
t3532 = t3415 * t3481 + t3421 * t3487;
t3336 = (-t3433 * t3532 - t3617) * t3486 - (g(3) * t3433 - t3479 * t3532) * t3480;
t3501 = 0.1e1 / qJ(3,3);
t3607 = t3336 * t3501;
t3531 = t3416 * t3483 + t3422 * t3489;
t3337 = (-t3434 * t3531 - t3617) * t3488 - (g(3) * t3434 - t3479 * t3531) * t3482;
t3502 = 0.1e1 / qJ(3,2);
t3606 = t3337 * t3502;
t3530 = t3417 * t3485 + t3423 * t3491;
t3338 = (-t3435 * t3530 - t3617) * t3490 - (g(3) * t3435 - t3479 * t3530) * t3484;
t3503 = 0.1e1 / qJ(3,1);
t3605 = t3338 * t3503;
t3387 = -t3436 * t3468 + t3442 * t3474;
t3604 = t3387 * t3467;
t3603 = t3387 * t3473;
t3388 = t3436 * t3474 + t3442 * t3468;
t3602 = t3388 * t3467;
t3601 = t3388 * t3473;
t3389 = -t3437 * t3470 + t3443 * t3476;
t3600 = t3389 * t3469;
t3599 = t3389 * t3475;
t3390 = t3437 * t3476 + t3443 * t3470;
t3598 = t3390 * t3469;
t3597 = t3390 * t3475;
t3391 = -t3438 * t3472 + t3444 * t3478;
t3596 = t3391 * t3471;
t3595 = t3391 * t3477;
t3392 = t3438 * t3478 + t3444 * t3472;
t3594 = t3392 * t3471;
t3593 = t3392 * t3477;
t3393 = -t3439 * t3481 + t3445 * t3487;
t3592 = t3393 * t3480;
t3591 = t3393 * t3486;
t3394 = t3439 * t3487 + t3445 * t3481;
t3590 = t3394 * t3480;
t3589 = t3394 * t3486;
t3395 = -t3440 * t3483 + t3446 * t3489;
t3588 = t3395 * t3482;
t3587 = t3395 * t3488;
t3396 = t3440 * t3489 + t3446 * t3483;
t3586 = t3396 * t3482;
t3585 = t3396 * t3488;
t3397 = -t3441 * t3485 + t3447 * t3491;
t3584 = t3397 * t3484;
t3583 = t3397 * t3490;
t3398 = t3441 * t3491 + t3447 * t3485;
t3582 = t3398 * t3484;
t3581 = t3398 * t3490;
t3562 = t3449 * t3450;
t3553 = t3387 * t3616;
t3552 = t3388 * t3616;
t3551 = t3389 * t3615;
t3550 = t3390 * t3615;
t3549 = t3391 * t3614;
t3548 = t3392 * t3614;
t3547 = t3393 * t3613;
t3546 = t3394 * t3613;
t3545 = t3395 * t3612;
t3544 = t3396 * t3612;
t3543 = t3397 * t3611;
t3542 = t3398 * t3611;
t3369 = t3449 * t3504 + t3452 * t3620;
t3541 = t3343 * t3387 + t3369 * t3388;
t3370 = t3449 * t3505 + t3452 * t3621;
t3540 = t3346 * t3389 + t3370 * t3390;
t3371 = t3449 * t3506 + t3452 * t3622;
t3539 = t3349 * t3391 + t3371 * t3392;
t3372 = t3449 * t3507 + t3452 * t3623;
t3538 = t3352 * t3393 + t3372 * t3394;
t3373 = t3449 * t3508 + t3452 * t3624;
t3537 = t3355 * t3395 + t3373 * t3396;
t3374 = t3449 * t3509 + t3452 * t3625;
t3536 = t3359 * t3397 + t3374 * t3398;
t3529 = t3448 * t3497 + t3451 * t3496;
t3368 = -t3490 * g(3) + t3484 * t3530;
t3367 = -g(3) * t3488 + t3482 * t3531;
t3366 = -g(3) * t3486 + t3480 * t3532;
t3365 = -g(3) * t3477 + t3471 * t3533;
t3364 = -t3475 * g(3) + t3469 * t3534;
t3363 = -g(3) * t3473 + t3467 * t3535;
t1 = [-g(1) * m(4) + (t3338 * t3583 - t3542) * t3503 + (t3337 * t3587 - t3544) * t3502 + (t3336 * t3591 - t3546) * t3501 + (t3335 * t3595 - t3548) * t3500 + (t3334 * t3599 - t3550) * t3499 + (t3333 * t3603 - t3552) * t3498 + (-t3363 * t3604 - t3364 * t3600 - t3365 * t3596 - t3366 * t3592 - t3367 * t3588 - t3368 * t3584) * m(3); -g(2) * m(4) + (t3338 * t3581 + t3543) * t3503 + (t3337 * t3585 + t3545) * t3502 + (t3336 * t3589 + t3547) * t3501 + (t3335 * t3593 + t3549) * t3500 + (t3334 * t3597 + t3551) * t3499 + (t3333 * t3601 + t3553) * t3498 + (-t3363 * t3602 - t3364 * t3598 - t3365 * t3594 - t3366 * t3590 - t3367 * t3586 - t3368 * t3582) * m(3); t3467 * t3610 + t3469 * t3609 + t3471 * t3608 + t3480 * t3607 + t3482 * t3606 + t3484 * t3605 - g(3) * m(4) + (t3363 * t3473 + t3364 * t3475 + t3365 * t3477 + t3366 * t3486 + t3367 * t3488 + t3368 * t3490) * m(3); -t3360 * t3503 * t3543 + (-t3359 * t3484 - t3360 * t3581) * t3605 - t3356 * t3502 * t3545 + (-t3355 * t3482 - t3356 * t3585) * t3606 - t3342 * t3501 * t3547 + (-t3342 * t3589 - t3352 * t3480) * t3607 - t3341 * t3500 * t3549 + (-t3341 * t3593 - t3349 * t3471) * t3608 - t3340 * t3499 * t3551 + (-t3340 * t3597 - t3346 * t3469) * t3609 - t3339 * t3498 * t3553 + (-t3339 * t3601 - t3343 * t3467) * t3610 - (-g(2) * t3619 + t3529 * g(3)) * t3453 + t3450 * (t3529 * g(2) + g(3) * t3619) + (-(t3359 * t3490 - t3360 * t3582) * t3368 - (t3355 * t3488 - t3356 * t3586) * t3367 - (-t3342 * t3590 + t3352 * t3486) * t3366 - (-t3341 * t3594 + t3349 * t3477) * t3365 - (-t3340 * t3598 + t3346 * t3475) * t3364 - (-t3339 * t3602 + t3343 * t3473) * t3363) * m(3); -t3626 * g(3) + (-t3360 * t3542 + (t3360 * t3583 - t3374 * t3484) * t3338) * t3503 + (-t3356 * t3544 + (t3356 * t3587 - t3373 * t3482) * t3337) * t3502 + (-t3342 * t3546 + (t3342 * t3591 - t3372 * t3480) * t3336) * t3501 + (-t3341 * t3548 + (t3341 * t3595 - t3371 * t3471) * t3335) * t3500 + (-t3340 * t3550 + (t3340 * t3599 - t3370 * t3469) * t3334) * t3499 + (-t3339 * t3552 + (t3339 * t3603 - t3369 * t3467) * t3333) * t3498 + (-t3529 * t3450 - t3453 * t3619) * g(1) + (-(t3360 * t3584 + t3374 * t3490) * t3368 - (t3356 * t3588 + t3373 * t3488) * t3367 - (t3342 * t3592 + t3372 * t3486) * t3366 - (t3341 * t3596 + t3371 * t3477) * t3365 - (t3340 * t3600 + t3370 * t3475) * t3364 - (t3339 * t3604 + t3369 * t3473) * t3363) * m(3); t3626 * g(2) + (-t3450 * t3554 + (t3453 * t3496 + t3497 * t3562) * t3451 + (t3453 * t3497 - t3496 * t3562) * t3448) * g(1) + (t3536 * t3338 * t3490 + (-t3359 * t3398 + t3374 * t3397) * t3611) * t3503 + (t3537 * t3337 * t3488 + (-t3355 * t3396 + t3373 * t3395) * t3612) * t3502 + (t3538 * t3336 * t3486 + (-t3352 * t3394 + t3372 * t3393) * t3613) * t3501 + (t3539 * t3335 * t3477 + (-t3349 * t3392 + t3371 * t3391) * t3614) * t3500 + (t3540 * t3334 * t3475 + (-t3346 * t3390 + t3370 * t3389) * t3615) * t3499 + (t3541 * t3333 * t3473 + (-t3343 * t3388 + t3369 * t3387) * t3616) * t3498 + (-t3363 * t3467 * t3541 - t3364 * t3469 * t3540 - t3365 * t3471 * t3539 - t3366 * t3480 * t3538 - t3367 * t3482 * t3537 - t3368 * t3484 * t3536) * m(3);];
taugX  = t1;
