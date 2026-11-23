// Optional owner diagnostics. Public model builds have no dependency on
// private implementations and incur no diagnostic work.
#pragma once
#define CART_TRACE(i, field, value) ((void)0)
#define CART_INITIAL(i, value) ((void)0)
